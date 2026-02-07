let running = false;

self.onmessage = (event) => {
  const message = event.data;
  if (message.type === "start") {
    running = true;
    const config = message.config;
    const trips = message.trips;
    const rng = createRng(message.seed || 1);
    optimize(config, trips, rng);
  }
  if (message.type === "stop") {
    running = false;
  }
};

function createRng(seed) {
  let state = seed >>> 0;
  return function rand() {
    state = (1664525 * state + 1013904223) >>> 0;
    return state / 4294967296;
  };
}

function randRange(rng, min, max) {
  return min + (max - min) * rng();
}

function randNormal(rng) {
  const u = Math.max(rng(), 1e-9);
  const v = Math.max(rng(), 1e-9);
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function distance(ax, ay, bx, by) {
  const dx = ax - bx;
  const dy = ay - by;
  return Math.hypot(dx, dy);
}

function segmentTime(length, vMax, accel, vStop) {
  const vmax = Math.max(vMax, vStop + 1e-6);
  const a = Math.max(accel, 1e-6);
  const vs = Math.min(vStop, vmax);
  const dAccel = (vmax * vmax - vs * vs) / (2 * a);
  if (2 * dAccel < length) {
    const tAccel = (vmax - vs) / a;
    const cruise = (length - 2 * dAccel) / vmax;
    return 2 * tAccel + cruise;
  }
  const vPeak = Math.sqrt(vs * vs + a * length);
  return 2 * (vPeak - vs) / a;
}

function initNetwork(config, rng) {
  const lineCount = Math.floor(randRange(rng, config.minLines, config.maxLines + 1));
  const weights = [];
  for (let i = 0; i < lineCount; i += 1) {
    weights.push(rng() + 0.2);
  }
  const weightSum = weights.reduce((a, b) => a + b, 0);
  const lines = weights.map((w) => ({
    length: (w / weightSum) * config.totalTrack,
    center: { x: 0, y: 0 },
    angle: randRange(rng, 0, Math.PI * 2),
    stations: [0, 1],
  }));

  lines.forEach((line) => {
    const half = line.length * 0.5;
    const dx = Math.abs(Math.cos(line.angle) * half);
    const dy = Math.abs(Math.sin(line.angle) * half);
    line.center.x = randRange(rng, dx, config.worldWidth - dx);
    line.center.y = randRange(rng, dy, config.worldHeight - dy);

    const interiorCount = Math.floor(randRange(rng, 0, Math.max(1, config.maxStations - 1)));
    for (let i = 0; i < interiorCount; i += 1) {
      line.stations.push(randRange(rng, 0.08, 0.92));
    }
    line.stations = normalizeStations(line.stations);
  });

  return {
    lines,
  };
}

function normalizeStations(stations) {
  const unique = new Set(stations.map((s) => clamp(s, 0, 1)));
  unique.add(0);
  unique.add(1);
  const sorted = Array.from(unique).sort((a, b) => a - b);
  if (sorted.length < 2) {
    return [0, 1];
  }
  return sorted;
}

function clampLine(line, config) {
  const half = line.length * 0.5;
  const dx = Math.abs(Math.cos(line.angle) * half);
  const dy = Math.abs(Math.sin(line.angle) * half);
  line.center.x = clamp(line.center.x, dx, config.worldWidth - dx);
  line.center.y = clamp(line.center.y, dy, config.worldHeight - dy);
}

function mutateNetwork(network, config, rng) {
  const next = {
    lines: network.lines.map((line) => ({
      length: line.length,
      center: { x: line.center.x, y: line.center.y },
      angle: line.angle,
      stations: [...line.stations],
    })),
  };

  const minLineLength = Math.max(40, config.totalTrack * 0.04);

  const roll = rng();
  if (roll < 0.2) {
    const line = next.lines[Math.floor(randRange(rng, 0, next.lines.length))];
    line.center.x += randNormal(rng) * (config.worldWidth * 0.03);
    line.center.y += randNormal(rng) * (config.worldHeight * 0.03);
    clampLine(line, config);
  } else if (roll < 0.35) {
    const line = next.lines[Math.floor(randRange(rng, 0, next.lines.length))];
    line.angle += randNormal(rng) * 0.2;
    clampLine(line, config);
  } else if (roll < 0.5 && next.lines.length > 1) {
    const i = Math.floor(randRange(rng, 0, next.lines.length));
    let j = Math.floor(randRange(rng, 0, next.lines.length));
    if (i === j) {
      j = (j + 1) % next.lines.length;
    }
    const lineA = next.lines[i];
    const lineB = next.lines[j];
    const delta = randRange(rng, -1, 1) * config.totalTrack * 0.05;
    if (lineA.length + delta > minLineLength && lineB.length - delta > minLineLength) {
      lineA.length += delta;
      lineB.length -= delta;
      clampLine(lineA, config);
      clampLine(lineB, config);
    }
  } else if (roll < 0.62) {
    mutateStations(next, config, rng);
  } else if (roll < 0.74) {
    if (next.lines.length < config.maxLines) {
      addLine(next, config, rng, minLineLength);
    }
  } else if (roll < 0.82) {
    if (next.lines.length > config.minLines) {
      removeLine(next, config, rng);
    }
  } else {
    const line = next.lines[Math.floor(randRange(rng, 0, next.lines.length))];
    line.center.x += randNormal(rng) * (config.worldWidth * 0.02);
    line.center.y += randNormal(rng) * (config.worldHeight * 0.02);
    line.angle += randNormal(rng) * 0.1;
    clampLine(line, config);
  }

  return next;
}

function addLine(network, config, rng, minLineLength) {
  if (!network.lines.length) {
    return;
  }
  const donorIndex = Math.floor(randRange(rng, 0, network.lines.length));
  const donor = network.lines[donorIndex];
  const take = Math.max(minLineLength, donor.length * 0.3);
  if (donor.length - take < minLineLength) {
    return;
  }
  donor.length -= take;

  const angle = randRange(rng, 0, Math.PI * 2);
  const half = take * 0.5;
  const dx = Math.abs(Math.cos(angle) * half);
  const dy = Math.abs(Math.sin(angle) * half);
  const center = {
    x: randRange(rng, dx, config.worldWidth - dx),
    y: randRange(rng, dy, config.worldHeight - dy),
  };

  const stations = normalizeStations([0, 1, randRange(rng, 0.2, 0.8)]);
  network.lines.push({ length: take, center, angle, stations });
}

function removeLine(network, config, rng) {
  if (network.lines.length <= 1) {
    return;
  }
  const idx = Math.floor(randRange(rng, 0, network.lines.length));
  const removed = network.lines.splice(idx, 1)[0];
  const receiver = network.lines[Math.floor(randRange(rng, 0, network.lines.length))];
  receiver.length += removed.length;
  clampLine(receiver, config);
}

function mutateStations(network, config, rng) {
  const line = network.lines[Math.floor(randRange(rng, 0, network.lines.length))];
  const interior = line.stations.slice(1, -1);
  const interiorCount = interior.length;

  const roll = rng();
  if (roll < 0.33 && line.stations.length < config.maxStations) {
    line.stations.push(randRange(rng, 0.08, 0.92));
  } else if (roll < 0.66 && interiorCount > 0) {
    const removeIndex = Math.floor(randRange(rng, 0, interiorCount));
    line.stations.splice(removeIndex + 1, 1);
  } else if (interiorCount > 0) {
    const moveIndex = Math.floor(randRange(rng, 0, interiorCount));
    const idx = moveIndex + 1;
    line.stations[idx] = clamp(line.stations[idx] + randNormal(rng) * 0.08, 0.05, 0.95);
  }

  line.stations = normalizeStations(line.stations);
}

function buildStations(network, config) {
  const stations = [];
  const stationLines = [];
  const lineStations = [];

  for (let lineIndex = 0; lineIndex < network.lines.length; lineIndex += 1) {
    const line = network.lines[lineIndex];
    const half = line.length * 0.5;
    const dx = Math.cos(line.angle) * half;
    const dy = Math.sin(line.angle) * half;
    const x1 = line.center.x - dx;
    const y1 = line.center.y - dy;
    const x2 = line.center.x + dx;
    const y2 = line.center.y + dy;

    const stationsOnLine = [];
    for (let s = 0; s < line.stations.length; s += 1) {
      const t = line.stations[s];
      const sx = x1 + (x2 - x1) * t;
      const sy = y1 + (y2 - y1) * t;

      let stationIndex = -1;
      for (let k = 0; k < stations.length; k += 1) {
        if (distance(sx, sy, stations[k].x, stations[k].y) <= config.mergeDistance) {
          stationIndex = k;
          break;
        }
      }
      if (stationIndex === -1) {
        stationIndex = stations.length;
        stations.push({ x: sx, y: sy });
        stationLines.push([]);
      }

      stationLines[stationIndex].push(lineIndex);
      stationsOnLine.push({ stationIndex, t, x: sx, y: sy });
    }

    stationsOnLine.sort((a, b) => a.t - b.t);
    lineStations.push(stationsOnLine);
  }

  return { stations, stationLines, lineStations };
}

function buildGraph(network, config) {
  const { stations, stationLines, lineStations } = buildStations(network, config);

  const nodes = [];
  const stationNodes = stations.map(() => []);

  for (let lineIndex = 0; lineIndex < lineStations.length; lineIndex += 1) {
    const entries = lineStations[lineIndex];
    for (let i = 0; i < entries.length; i += 1) {
      const entry = entries[i];
      const nodeId = nodes.length;
      nodes.push({ lineIndex, stationIndex: entry.stationIndex, x: entry.x, y: entry.y });
      stationNodes[entry.stationIndex].push(nodeId);
      entry.nodeId = nodeId;
    }
  }

  const adjacency = nodes.map(() => []);
  for (let lineIndex = 0; lineIndex < lineStations.length; lineIndex += 1) {
    const entries = lineStations[lineIndex];
    for (let i = 0; i < entries.length - 1; i += 1) {
      const a = entries[i];
      const b = entries[i + 1];
      const length = distance(a.x, a.y, b.x, b.y);
      const time = segmentTime(length, config.trainSpeed, config.accel, config.stopSpeed);
      adjacency[a.nodeId].push({ to: b.nodeId, weight: time });
      adjacency[b.nodeId].push({ to: a.nodeId, weight: time });
    }
  }

  for (let s = 0; s < stationNodes.length; s += 1) {
    const nodesAtStation = stationNodes[s];
    if (nodesAtStation.length > 1) {
      for (let i = 0; i < nodesAtStation.length; i += 1) {
        for (let j = i + 1; j < nodesAtStation.length; j += 1) {
          const a = nodesAtStation[i];
          const b = nodesAtStation[j];
          adjacency[a].push({ to: b, weight: config.transferPenalty });
          adjacency[b].push({ to: a, weight: config.transferPenalty });
        }
      }
    }
  }

  const distToStation = nodes.map(() => new Array(stations.length).fill(Infinity));
  for (let i = 0; i < nodes.length; i += 1) {
    const dist = dijkstra(i, adjacency);
    for (let s = 0; s < stations.length; s += 1) {
      let best = Infinity;
      const nodesAt = stationNodes[s];
      for (let k = 0; k < nodesAt.length; k += 1) {
        const d = dist[nodesAt[k]];
        if (d < best) {
          best = d;
        }
      }
      distToStation[i][s] = best;
    }
  }

  return { stations, stationNodes, distToStation };
}

function dijkstra(start, adjacency) {
  const n = adjacency.length;
  const dist = new Array(n).fill(Infinity);
  const visited = new Array(n).fill(false);
  dist[start] = 0;

  for (let i = 0; i < n; i += 1) {
    let best = Infinity;
    let node = -1;
    for (let j = 0; j < n; j += 1) {
      if (!visited[j] && dist[j] < best) {
        best = dist[j];
        node = j;
      }
    }
    if (node === -1) {
      break;
    }
    visited[node] = true;
    const edges = adjacency[node];
    for (let e = 0; e < edges.length; e += 1) {
      const edge = edges[e];
      const nextDist = dist[node] + edge.weight;
      if (nextDist < dist[edge.to]) {
        dist[edge.to] = nextDist;
      }
    }
  }

  return dist;
}

function evaluate(network, config, trips) {
  const graph = buildGraph(network, config);
  const { stations, stationNodes, distToStation } = graph;
  let total = 0;

  for (let i = 0; i < trips.length; i += 1) {
    const trip = trips[i];
    const direct = distance(trip.ax, trip.ay, trip.bx, trip.by) / config.walkSpeed;
    let best = direct;

    const walkToDest = stations.map((s) => distance(trip.bx, trip.by, s.x, s.y) / config.walkSpeed);

    for (let s = 0; s < stations.length; s += 1) {
      const walkStart = distance(trip.ax, trip.ay, stations[s].x, stations[s].y) / config.walkSpeed;
      if (walkStart + config.boardPenalty >= best) {
        continue;
      }
      const nodesAt = stationNodes[s];
      for (let n = 0; n < nodesAt.length; n += 1) {
        const nodeId = nodesAt[n];
        for (let t = 0; t < stations.length; t += 1) {
          const route = distToStation[nodeId][t];
          if (!Number.isFinite(route)) {
            continue;
          }
          const candidate = walkStart + config.boardPenalty + route + walkToDest[t];
          if (candidate < best) {
            best = candidate;
          }
        }
      }
    }

    total += best;
  }

  return total / trips.length;
}

function optimize(config, trips, rng) {
  let current = initNetwork(config, rng);
  let currentScore = evaluate(current, config, trips);
  let best = current;
  let bestScore = currentScore;

  const startTemp = 2.0;
  const endTemp = 0.15;
  const iterations = config.iterations;
  const reportEvery = Math.max(20, Math.floor(iterations / 150));

  for (let iter = 0; iter < iterations && running; iter += 1) {
    const temp = startTemp * Math.pow(endTemp / startTemp, iter / iterations);
    const candidate = mutateNetwork(current, config, rng);
    const score = evaluate(candidate, config, trips);
    const delta = score - currentScore;
    if (delta < 0 || Math.exp(-delta / temp) > rng()) {
      current = candidate;
      currentScore = score;
      if (score < bestScore) {
        best = candidate;
        bestScore = score;
      }
    }

    if (iter % reportEvery === 0 || iter === iterations - 1) {
      const bestStationsCount = best.lines.reduce((sum, line) => sum + line.stations.length, 0);
      self.postMessage({
        type: "update",
        iter,
        currentScore,
        bestScore,
        bestLines: best.lines.length,
        bestStations: bestStationsCount,
        bestNetwork: best,
      });
    }
  }

  self.postMessage({ type: "done" });
}
