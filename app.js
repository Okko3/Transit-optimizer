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
  if (!trips.length || !lastConfig) {
    return;
  }
  ctx.save();
  ctx.translate(view.offsetX, view.offsetY);
  ctx.scale(view.scale, view.scale);
  ctx.strokeStyle = "rgba(143,179,255,0.22)";
  ctx.lineWidth = 1 / view.scale;
  const sample = Math.min(trips.length, 60);
  for (let i = 0; i < sample; i += 1) {
    const t = trips[i];
    ctx.beginPath();
    ctx.moveTo(t.ax, t.ay);
    ctx.lineTo(t.bx, t.by);
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
