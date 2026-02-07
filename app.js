const form = document.getElementById("configForm");
const statusText = document.getElementById("statusText");
const bestTimeEl = document.getElementById("bestTime");
const currentTimeEl = document.getElementById("currentTime");
const bestLinesEl = document.getElementById("bestLines");
const bestStationsEl = document.getElementById("bestStations");
const iterTextEl = document.getElementById("iterText");
const startBtn = document.getElementById("startBtn");
const stopBtn = document.getElementById("stopBtn");
const resampleBtn = document.getElementById("resampleBtn");
const canvas = document.getElementById("map");
const ctx = canvas.getContext("2d");

let workers = [];
let bestNetwork = null;
let bestScore = Infinity;
let lastConfig = null;
let trips = [];
let tripPaths = [];
let runState = "Idle";
let activeWorkers = 0;
let view = {
  scale: 1,
  offsetX: 0,
  offsetY: 0,
  worldWidth: 5000,
  worldHeight: 5000,
  zoom: 1,
};

const palette = [
  "#f6c177",
  "#4cc9f0",
  "#f28482",
  "#90be6d",
  "#b5179e",
  "#ffd166",
  "#70d6ff",
  "#cdb4db",
];
const TRIP_DRAW_LIMIT = 60;

function readNumber(id, fallback) {
  const el = document.getElementById(id);
  const value = parseFloat(el.value);
  return Number.isFinite(value) ? value : fallback;
}

function createRng(seed) {
  let state = seed >>> 0;
  return function rand() {
    state = (1664525 * state + 1013904223) >>> 0;
    return state / 4294967296;
  };
}

function sampleTrips(config) {
  const seed = Math.floor(config.seed || 1);
  const rng = createRng(seed);
  const list = [];
  for (let i = 0; i < config.tripCount; i += 1) {
    const ax = rng() * config.worldWidth;
    const ay = rng() * config.worldHeight;
    const bx = rng() * config.worldWidth;
    const by = rng() * config.worldHeight;
    list.push({ ax, ay, bx, by });
  }
  return list;
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

  return { stations, stationNodes, nodes, adjacency };
}

function dijkstraWithPrev(start, adjacency) {
  const n = adjacency.length;
  const dist = new Array(n).fill(Infinity);
  const visited = new Array(n).fill(false);
  const prev = new Array(n).fill(-1);
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
        prev[edge.to] = node;
      }
    }
  }

  return { dist, prev };
}

function reconstructPath(prev, start, target) {
  const path = [];
  let cur = target;
  while (cur !== -1) {
    path.push(cur);
    if (cur === start) {
      break;
    }
    cur = prev[cur];
  }
  if (path[path.length - 1] !== start) {
    return null;
  }
  path.reverse();
  return path;
}

function buildTripPaths(tripList, network, config) {
  const count = Math.min(tripList.length, TRIP_DRAW_LIMIT);
  if (!network || !network.lines || !network.lines.length) {
    return tripList.slice(0, count).map((trip) => ({
      usesTrain: false,
      path: [
        { x: trip.ax, y: trip.ay },
        { x: trip.bx, y: trip.by },
      ],
    }));
  }

  const graph = buildGraph(network, config);
  const { stations, stationNodes, nodes, adjacency } = graph;
  if (!stations.length || !nodes.length) {
    return tripList.slice(0, count).map((trip) => ({
      usesTrain: false,
      path: [
        { x: trip.ax, y: trip.ay },
        { x: trip.bx, y: trip.by },
      ],
    }));
  }

  const allDijkstra = nodes.map((_, idx) => dijkstraWithPrev(idx, adjacency));
  const results = [];

  for (let i = 0; i < count; i += 1) {
    const trip = tripList[i];
    const directTime = distance(trip.ax, trip.ay, trip.bx, trip.by) / config.walkSpeed;
    let bestTime = directTime;
    let bestEntry = -1;
    let bestExit = -1;
    let bestPrev = null;

    const walkEnd = stations.map((s) => distance(trip.bx, trip.by, s.x, s.y) / config.walkSpeed);

    for (let s = 0; s < stations.length; s += 1) {
      const walkStart = distance(trip.ax, trip.ay, stations[s].x, stations[s].y) / config.walkSpeed;
      const walkStartTotal = walkStart + config.boardPenalty;
      if (walkStartTotal >= bestTime) {
        continue;
      }
      const entryNodes = stationNodes[s];
      for (let e = 0; e < entryNodes.length; e += 1) {
        const entryNode = entryNodes[e];
        const { dist, prev } = allDijkstra[entryNode];
        for (let exitNode = 0; exitNode < nodes.length; exitNode += 1) {
          const trainDist = dist[exitNode];
          if (!Number.isFinite(trainDist)) {
            continue;
          }
          const exitStation = nodes[exitNode].stationIndex;
          const total = walkStartTotal + trainDist + walkEnd[exitStation];
          if (total < bestTime) {
            bestTime = total;
            bestEntry = entryNode;
            bestExit = exitNode;
            bestPrev = prev;
          }
        }
      }
    }

    if (bestEntry >= 0 && bestExit >= 0 && bestPrev) {
      const trainPath = reconstructPath(bestPrev, bestEntry, bestExit);
      if (trainPath && trainPath.length) {
        const trainPoints = trainPath.map((idx) => ({ x: nodes[idx].x, y: nodes[idx].y }));
        results.push({
          usesTrain: true,
          path: [{ x: trip.ax, y: trip.ay }, ...trainPoints, { x: trip.bx, y: trip.by }],
        });
        continue;
      }
    }

    results.push({
      usesTrain: false,
      path: [
        { x: trip.ax, y: trip.ay },
        { x: trip.bx, y: trip.by },
      ],
    });
  }

  return results;
}

function readConfig() {
  const worldWidth = readNumber("worldWidth", 5000);
  const worldHeight = readNumber("worldHeight", 5000);
  const minLines = Math.max(1, Math.floor(readNumber("minLines", 1)));
  const maxLines = Math.max(minLines, Math.floor(readNumber("maxLines", 10)));
  const totalTrack = Math.max(50, readNumber("totalTrack", 10000));

  return {
    worldWidth,
    worldHeight,
    totalTrack,
    minLines,
    maxLines,
    trainSpeed: Math.max(0.1, readNumber("trainSpeed", 70)),
    walkSpeed: Math.max(0.1, readNumber("walkSpeed", 5)),
    accel: Math.max(0.01, readNumber("accel", 1)),
    stopSpeed: Math.max(0, readNumber("stopSpeed", 0)),
    boardPenalty: Math.max(0, readNumber("boardPenalty", 2)),
    transferPenalty: Math.max(0, readNumber("transferPenalty", 2)),
    tripCount: Math.max(10, Math.floor(readNumber("tripCount", 1000))),
    iterations: Math.max(200, Math.floor(readNumber("iterations", 30000))),
    maxStations: Math.max(2, Math.floor(readNumber("maxStations", 8))),
    mergeDistance: Math.max(0, readNumber("mergeDistance", 8)),
    workers: Math.max(0, Math.floor(readNumber("workers", 12))),
    seed: Math.floor(readNumber("seed", 1234)),
  };
}

function resizeCanvas() {
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.scale(dpr, dpr);
  updateView();
  draw();
}

function updateView() {
  if (!lastConfig) {
    return;
  }
  view.worldWidth = lastConfig.worldWidth;
  view.worldHeight = lastConfig.worldHeight;
  const rect = canvas.getBoundingClientRect();
  const scaleX = rect.width / view.worldWidth;
  const scaleY = rect.height / view.worldHeight;
  view.scale = Math.min(scaleX, scaleY) * view.zoom;
  view.offsetX = (rect.width - view.worldWidth * view.scale) * 0.5;
  view.offsetY = (rect.height - view.worldHeight * view.scale) * 0.5;
}

function screenToWorld(x, y) {
  return {
    x: (x - view.offsetX) / view.scale,
    y: (y - view.offsetY) / view.scale,
  };
}

function drawGrid() {
  const rect = canvas.getBoundingClientRect();
  ctx.fillStyle = "#0b1118";
  ctx.fillRect(0, 0, rect.width, rect.height);
  if (!lastConfig) {
    return;
  }
  ctx.save();
  ctx.translate(view.offsetX, view.offsetY);
  ctx.scale(view.scale, view.scale);
  ctx.strokeStyle = "rgba(255,255,255,0.05)";
  ctx.lineWidth = 1 / view.scale;
  const step = 100;
  for (let x = 0; x <= view.worldWidth; x += step) {
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, view.worldHeight);
    ctx.stroke();
  }
  for (let y = 0; y <= view.worldHeight; y += step) {
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(view.worldWidth, y);
    ctx.stroke();
  }
  ctx.strokeStyle = "rgba(246,193,119,0.18)";
  ctx.lineWidth = 2 / view.scale;
  ctx.strokeRect(0, 0, view.worldWidth, view.worldHeight);
  ctx.restore();
}

function drawTrips() {
  if (!tripPaths.length || !lastConfig) {
    return;
  }
  ctx.save();
  ctx.translate(view.offsetX, view.offsetY);
  ctx.scale(view.scale, view.scale);
  ctx.strokeStyle = "rgba(143,179,255,0.22)";
  ctx.lineWidth = 1 / view.scale;
  for (let i = 0; i < tripPaths.length; i += 1) {
    const t = tripPaths[i];
    const path = t.path;
    if (!path || path.length < 2) {
      continue;
    }
    ctx.beginPath();
    ctx.moveTo(path[0].x, path[0].y);
    for (let p = 1; p < path.length; p += 1) {
      ctx.lineTo(path[p].x, path[p].y);
    }
    ctx.stroke();
  }
  ctx.restore();
}

function drawNetwork(network) {
  if (!network) {
    return;
  }
  const lines = network.lines;
  ctx.save();
  ctx.translate(view.offsetX, view.offsetY);
  ctx.scale(view.scale, view.scale);
  lines.forEach((line, idx) => {
    const color = palette[idx % palette.length];
    const dx = Math.cos(line.angle) * (line.length * 0.5);
    const dy = Math.sin(line.angle) * (line.length * 0.5);
    const x1 = line.center.x - dx;
    const y1 = line.center.y - dy;
    const x2 = line.center.x + dx;
    const y2 = line.center.y + dy;

    ctx.strokeStyle = color;
    ctx.lineWidth = 3 / view.scale;
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();

    ctx.fillStyle = color;
    line.stations.forEach((t) => {
      const sx = x1 + (x2 - x1) * t;
      const sy = y1 + (y2 - y1) * t;
      ctx.beginPath();
      ctx.arc(sx, sy, 4 / view.scale, 0, Math.PI * 2);
      ctx.fill();
    });
  });
  ctx.restore();
}

function draw() {
  drawGrid();
  drawTrips();
  drawNetwork(bestNetwork);
}

function updateStats(stats) {
  if (stats.bestScore != null) {
    bestTimeEl.textContent = stats.bestScore.toFixed(2);
  }
  if (stats.currentScore != null) {
    currentTimeEl.textContent = stats.currentScore.toFixed(2);
  }
  if (stats.bestLines != null) {
    bestLinesEl.textContent = stats.bestLines.toString();
  }
  if (stats.bestStations != null) {
    bestStationsEl.textContent = stats.bestStations.toString();
  }
  if (stats.iter != null) {
    iterTextEl.textContent = stats.iter.toString();
  }
}

function clearWorkers() {
  workers.forEach((worker) => worker.terminate());
  workers = [];
}

function startOptimization() {
  stopOptimization();
  const config = readConfig();
  lastConfig = config;
  trips = sampleTrips(config);
  tripPaths = buildTripPaths(trips, bestNetwork, config);
  bestScore = Infinity;
  bestNetwork = null;
  updateStats({ bestScore: null, currentScore: null, bestLines: null, bestStations: null, iter: null });
  draw();

  const workerCount = config.workers > 0
    ? config.workers
    : Math.max(1, Math.min(navigator.hardwareConcurrency || 4, 12));

  runState = `Running with ${workerCount} worker${workerCount === 1 ? "" : "s"}`;
  statusText.textContent = runState;
  startBtn.disabled = true;
  stopBtn.disabled = false;
  activeWorkers = workerCount;

  for (let i = 0; i < workerCount; i += 1) {
    const worker = new Worker("worker.js");
    worker.onmessage = (event) => handleWorkerMessage(event.data);
    worker.postMessage({
      type: "start",
      config,
      trips,
      seed: config.seed + i * 1013,
    });
    workers.push(worker);
  }
}

function stopOptimization() {
  clearWorkers();
  startBtn.disabled = false;
  stopBtn.disabled = true;
  runState = "Idle";
  statusText.textContent = runState;
  activeWorkers = 0;
}

function handleWorkerMessage(message) {
  if (message.type === "update") {
    if (message.bestScore < bestScore) {
      bestScore = message.bestScore;
      bestNetwork = message.bestNetwork;
    }
    if (message.trips) {
      trips = message.trips;
    }
    if (lastConfig) {
      tripPaths = buildTripPaths(trips, bestNetwork, lastConfig);
    }
    updateStats({
      bestScore,
      currentScore: message.currentScore,
      bestLines: message.bestLines,
      bestStations: message.bestStations,
      iter: message.iter,
    });
    draw();
  }
  if (message.type === "done") {
    activeWorkers = Math.max(0, activeWorkers - 1);
    if (activeWorkers === 0) {
      runState = "Complete";
      statusText.textContent = runState;
    }
  }
}

startBtn.addEventListener("click", () => {
  startOptimization();
});

stopBtn.addEventListener("click", () => {
  stopOptimization();
});

resampleBtn.addEventListener("click", () => {
  lastConfig = readConfig();
  trips = sampleTrips(lastConfig);
  tripPaths = buildTripPaths(trips, bestNetwork, lastConfig);
  draw();
});

window.addEventListener("resize", resizeCanvas);

canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  const delta = Math.sign(event.deltaY) * -0.05;
  view.zoom = Math.min(2.5, Math.max(0.6, view.zoom + delta));
  updateView();
  draw();
});

canvas.addEventListener("mousemove", (event) => {
  if (!lastConfig) {
    return;
  }
  const rect = canvas.getBoundingClientRect();
  const pos = screenToWorld(event.clientX - rect.left, event.clientY - rect.top);
  statusText.textContent = `Hover: ${pos.x.toFixed(0)}, ${pos.y.toFixed(0)}`;
});

canvas.addEventListener("mouseleave", () => {
  statusText.textContent = runState;
});

lastConfig = readConfig();
trips = sampleTrips(lastConfig);
resizeCanvas();
