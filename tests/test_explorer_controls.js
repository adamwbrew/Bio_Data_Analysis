/* Exercise the actual offline explorer script with synthetic measured grids.
 * This is a small DOM/Plotly contract harness, not browser automation or a
 * second implementation of the filtering calculation.
 */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

class Element {
  constructor(tag = 'div') {
    this.tagName = tag.toUpperCase(); this.value = ''; this.textContent = '';
    this.innerHTML = ''; this.hidden = false; this.disabled = false;
    this.checked = false; this.children = []; this.dataset = {}; this.attributes = {}; this.events = new Map();
    this.style = {setProperty() {}};
    this.classList = {toggle() {}, add() {}, remove() {}};
  }
  get options() { return this.children; }
  appendChild(child) { this.children.push(child); return child; }
  append(...children) { this.children.push(...children); }
  add(option) { this.children.push(option); if (this.children.length === 1) this.value = option.value; }
  replaceChildren(...children) { this.children = [...children]; if (this.tagName === 'SELECT') this.value = ''; }
  setAttribute(name, value) { this.attributes[name] = String(value); }
  addEventListener(name, handler) { this['on' + name] = handler; }
  on(name, handler) { const handlers = this.events.get(name) || []; handlers.push(handler); this.events.set(name, handlers); return this; }
  emit(name, event) { (this.events.get(name) || []).forEach(handler => handler(event)); }
}

function grid(width) {
  return {width, x: [width / 2, 1.5 * width, 2.5 * width], y: [width / 2, 1.5 * width],
    n_pla: [[0, 4, 5], [10, 5, 1]], n_chat: [[2, 2, 2], [2, 0, 2]],
    n_all: [[2, 6, 7], [12, 5, 3]], n_positive: [[0, 2, 0], [5, 2, 1]],
    pla_area: [[0, 8, 15], [20, 10, 1]], n_large: [[0, 1, 0], [4, 1, 0]],
    overlap_area: [[0, 1, 0], [2, .5, .2]]};
}

function harness() {
  const template = fs.readFileSync(path.join(__dirname, '../bio_visuals/explorer.html'), 'utf8');
  const elements = new Map(), tabs = [], root = new Element();
  for (const match of template.matchAll(/<(\w+)\b([^>]*\bdata-role="([^"]+)"[^>]*)>/g)) {
    const [, tag, attributes, role] = match, node = new Element(tag);
    node.value = (attributes.match(/\bvalue="([^"]*)"/) || [])[1] || '';
    node.checked = /\bchecked\b/.test(attributes);
    node.hidden = /\bhidden\b/.test(attributes);
    if (tag === 'select') {
      const contents = template.slice(match.index + match[0].length).split('</select>')[0];
      for (const option of contents.matchAll(/<option([^>]*)>([^<]*)<\/option>/g)) {
        const child = new Element('option'); child.value = (option[1].match(/\bvalue="([^"]*)"/) || [])[1] || option[2];
        child.textContent = option[2]; node.children.push(child);
        if (node.children.length === 1 || /\bselected\b/.test(option[1])) node.value = child.value;
      }
    }
    elements.set(role, node);
  }
  for (const match of template.matchAll(/<button[^>]*data-view="([^"]+)"[^>]*>/g)) {
    const node = new Element('button'); node.dataset.view = match[1]; tabs.push(node);
  }
  root.querySelector = selector => {
    const role = (selector.match(/data-role=["']([^"']+)["']/) || [])[1];
    assert.ok(elements.has(role), 'Harness must represent actual control ' + selector);
    return elements.get(role);
  };
  root.querySelectorAll = selector => { assert.equal(selector, '[data-view]'); return tabs; };
  const payload = {sections: ['21F', '22F'].map((animal, order) => ({
    job: String(order + 1), animal, order, group: order ? 'Running' : 'Home cage',
    arm: 'TDP43–HK1', color: order ? '#d06738' : '#087e8b', ready: false,
    candidate_center: null, proposed: null, fine: grid(100),
    grids: {'250': grid(250), '500': grid(500), '1000': grid(1000)}
  })), colorscales: {Magma: [[0, '#000'], [1, '#fff']], Viridis: [[0, '#000'], [1, '#fff']], Cividis: [[0, '#000'], [1, '#fff']]}};
  let timerId = 0; const timers = new Map(), plots = [];
  const context = {
    document: {getElementById: () => root, createElement: tag => new Element(tag), createTextNode: text => ({textContent: text})},
    Option: function(text, value) { this.textContent = text; this.value = value; },
    Plotly: {react: async (node, traces, layout) => { const rendered = JSON.parse(JSON.stringify({traces, layout})); plots.push(rendered); node.layout = rendered.layout; node._fullLayout = rendered.layout; return node; }},
    setTimeout: callback => { timers.set(++timerId, callback); return timerId; },
    clearTimeout: id => timers.delete(id), console
  };
  const scripts = [...template.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(match => match[1]);
  const script = scripts.find(text => text.includes('const root=document.getElementById'));
  assert.ok(script, 'Must test the actual explorer inline script');
  vm.runInNewContext(script.replace('__DATA__', JSON.stringify(payload)), context, {filename: 'explorer.html'});
  return {elements, tabs, plots, timers,
    async settle() { for (let i = 0; i < 12; i++) await Promise.resolve(); },
    flushTimers() { const pending = [...timers.values()]; timers.clear(); pending.forEach(callback => callback()); },
    change(role, value) { const node = elements.get(role); node.value = String(value); assert.equal(typeof node.onchange, 'function', role + ' change handler'); node.onchange(); },
    input(role, value) { const node = elements.get(role); node.value = String(value); assert.equal(typeof node.oninput, 'function', role + ' input handler'); node.oninput(); },
    click(role) { const node = elements.get(role); assert.ok(node, role + ' control exists'); assert.equal(typeof node.onclick, 'function', role + ' click handler'); node.onclick(); },
    relayout(event) { assert.ok(elements.get('plot').events.has('plotly_relayout'), 'Manual camera changes are observed'); elements.get('plot').emit('plotly_relayout', event); },
    last() { assert.ok(plots.length); return plots[plots.length - 1]; }
  };
}

const near = (actual, expected, message) => assert.ok(Math.abs(actual - expected) < 1e-9, message + ': ' + actual + ' versus ' + expected);
const relativeEye = camera => Object.fromEntries(['x', 'y', 'z'].map(axis => [axis, camera.eye[axis] - ((camera.center || {})[axis] || 0)]));
const distance = camera => { const e = relativeEye(camera); return Math.hypot(e.x, e.y, e.z); };
const elevation = camera => { const e = relativeEye(camera); return Math.atan2(e.z, Math.hypot(e.x, e.y)); };
function finiteCamera(camera) {
  for (const vector of [camera.eye, camera.up, camera.center]) {
    if (vector) for (const axis of ['x', 'y', 'z']) assert.ok(Number.isFinite(vector[axis]), 'Finite camera ' + axis);
  }
  const e = relativeEye(camera), u = camera.up;
  if (u) assert.ok(Math.hypot(e.y*u.z-e.z*u.y, e.z*u.x-e.x*u.z, e.x*u.y-e.y*u.x) > 1e-9, 'Camera up must not be parallel to viewing axis');
}

async function cameraChecks(app) {
  const get = name => app.elements.get(name), camera = () => app.last().layout.scene.camera;
  assert.equal(get('camera-controls').hidden, false, '3D measurement view exposes camera controls');
  const measuredTraces = JSON.stringify(app.last().traces);
  const defaultEye = {...camera().eye};
  for (const [role, axis, sign] of [
    ['camera-top', 'z', 1], ['camera-bottom', 'z', -1],
    ['camera-x-positive', 'x', 1], ['camera-x-negative', 'x', -1],
    ['camera-y-positive', 'y', 1], ['camera-y-negative', 'y', -1]
  ]) {
    app.click(role); await app.settle(); const c = camera(); finiteCamera(c);
    if (role === 'camera-top') assert.equal(app.last().layout.scene.zaxis.visible, false, 'Exact top view hides the edge-on z labels');
    assert.equal(c.projection.type, 'orthographic', role + ' gives a parallel axis view');
    const range = app.last().layout.scene[axis + 'axis'].range;
    const axisDirection = range && range[1] < range[0] ? -1 : 1;
    assert.ok(c.eye[axis] * sign * axisDirection > 0, role + ' looks from the named data-coordinate side, including reversed image y');
    for (const other of ['x', 'y', 'z'].filter(a => a !== axis)) near(c.eye[other], 0, role + ' is axis aligned');
  }
  app.click('camera-oblique'); await app.settle();
  assert.equal(camera().projection.type, 'perspective', 'Oblique restores perspective');
  for (const axis of ['x', 'y', 'z']) assert.notEqual(camera().eye[axis], 0, 'Oblique shows all three dimensions');

  const manuallyDragged = {eye: {x: 2, y: 1, z: .5}, center: {x: .1, y: -.2, z: .3},
    up: {x: 0, y: 0, z: 1}, projection: {type: 'orthographic'}};
  app.click('camera-top'); await app.settle();
  app.relayout({'scene.camera': manuallyDragged});
  await app.settle();
  for (const axis of ['x', 'y', 'z']) assert.equal(app.last().layout.scene[axis + 'axis'].visible, true, 'Manual movement away from a fixed axis view restores ' + axis + ' labels');
  app.click('rotate-left'); await app.settle(); const turned = camera();
  near(distance(turned), distance(manuallyDragged), 'Turn keeps manually selected zoom distance');
  near(turned.eye.z, manuallyDragged.eye.z, 'Azimuth turn preserves vertical eye position');
  const afterEye = relativeEye(turned), beforeEye = relativeEye(manuallyDragged);
  const delta = Math.atan2(afterEye.y, afterEye.x) - Math.atan2(beforeEye.y, beforeEye.x);
  near(Math.abs(delta), Math.PI / 6, 'Turn advances 30 degrees from manually dragged camera');
  assert.deepEqual(turned.center, manuallyDragged.center, 'Turn preserves the manually panned center');
  assert.deepEqual(turned.projection, manuallyDragged.projection, 'Turn preserves projection mode');
  app.click('rotate-right'); await app.settle();
  for (const axis of ['x', 'y', 'z']) near(camera().eye[axis], manuallyDragged.eye[axis], 'Inverse turn restores ' + axis);
  app.click('tilt-up'); await app.settle();
  near(elevation(camera()) - elevation(manuallyDragged), Math.PI / 12, 'Tilt advances 15 degrees');
  near(distance(camera()), distance(manuallyDragged), 'Tilt preserves zoom distance');
  assert.deepEqual(camera().center, manuallyDragged.center, 'Tilt preserves center');
  app.click('tilt-down'); await app.settle();

  app.relayout({'scene.camera.eye.x': 1.5, 'scene.camera.center.y': .25});
  app.click('rotate-left'); await app.settle();
  near(distance(camera()), Math.hypot(1.5 - .1, 1 - .25, .5 - .3), 'Dotted relayout updates merge with previous manual camera');
  near(camera().center.y, .25, 'Dotted center update survives a turn');
  near(camera().center.x, .1, 'Dotted update keeps unchanged center coordinates');

  app.click('camera-top'); await app.settle(); const topDistance = distance(camera()), topUp = {...camera().up};
  app.click('rotate-left'); await app.settle(); finiteCamera(camera());
  assert.notDeepEqual(camera().up, topUp, 'Turning an exact top view rotates the view instead of becoming a no-op');
  near(elevation(camera()), Math.PI / 2, 'Turning at a pole keeps the top axis view');
  app.click('tilt-down'); await app.settle(); finiteCamera(camera());
  near(distance(camera()), topDistance, 'Tilt away from a pole preserves camera distance');
  assert.ok(Math.hypot(camera().eye.x, camera().eye.y) > 0, 'Top view can tilt away from the pole');
  for (let i = 0; i < 14; i++) { app.click('tilt-down'); await app.settle(); }
  finiteCamera(camera()); near(elevation(camera()), -Math.PI / 2, 'Tilt clamps at the bottom pole');
  app.click('camera-reset'); await app.settle();
  for (const axis of ['x', 'y', 'z']) near(camera().eye[axis], defaultEye[axis], 'Reset restores default eye ' + axis);
  assert.equal(JSON.stringify(app.last().traces), measuredTraces, 'Camera changes never alter measured coordinates or values');
  assert.equal(get('plot').events.get('plotly_relayout').length, 1, 'Repeated renders do not duplicate manual-camera listeners');
}

(async function run() {
  const app = harness(), get = name => app.elements.get(name);
  await app.settle();
  assert.equal(get('minimum').disabled, true, 'Counts do not use a PLA denominator cutoff');
  assert.ok(get('minimum-help').textContent.trim(), 'Inactive slider needs an adjacent explanation');
  assert.equal(get('minimum-status').hidden, true, 'Counts do not display an applied fraction cutoff');
  await cameraChecks(app);

  app.change('metric', 'fraction'); await app.settle();
  assert.equal(get('minimum').disabled, false);
  const beforeInput = app.plots.length;
  app.input('minimum', 5);
  assert.match(get('minimum-value').textContent, /5/, 'Readout changes during drag, before plot redraw');
  assert.equal(app.plots.length, beforeInput, 'Plot updates are throttled');
  app.flushTimers(); await app.settle();
  const expectedFraction = [[null, null, 0], [50, 40, null]];
  assert.deepEqual(app.last().traces[0].z, expectedFraction,
    'Fraction masks below cutoff, includes exact cutoff, preserves supported zero and excludes zero-PLA bins');
  assert.equal(get('minimum-status').hidden, false);
  assert.match(get('minimum-status').textContent, /3 \/ 6/, 'Applied cutoff reports supported versus observed bins');

  app.input('minimum', 1); app.input('minimum', 10);
  assert.equal(app.timers.size, 1, 'Fast drag changes coalesce to the latest cutoff');
  app.change('minimum', 10); await app.settle();
  assert.equal(app.timers.size, 0, 'Change flushes the pending drag update');
  assert.deepEqual(app.last().traces[0].z, [[null, null, null], [50, null, null]]);

  app.change('metric', 'mean_area'); app.change('minimum', 5); await app.settle();
  assert.deepEqual(app.last().traces[0].z, [[null, null, 3], [2, 2, null]], 'Mean area uses the same PLA denominator threshold');

  app.change('metric', 'n_pla'); await app.settle();
  assert.equal(get('minimum').disabled, true);
  assert.deepEqual(app.last().traces[0].z, grid(500).n_pla.map(row => row.map(Math.log1p)),
    'Raw counts retain observed bins including supported zero PLA, regardless of cutoff');

  app.tabs.find(tab => tab.dataset.view === 'layers').onclick();
  app.change('layer-view', 'overlay'); app.change('layer-metric', 'fraction');
  app.input('minimum', 5); app.flushTimers(); await app.settle();
  assert.equal(get('minimum').disabled, false);
  assert.equal(get('camera-controls').hidden, true, 'Flat layers hide 3D camera controls');
  const layerMaps = app.last().traces.filter(trace => trace.type === 'heatmap');
  assert.equal(layerMaps.length, 2);
  layerMaps.forEach(trace => assert.deepEqual(trace.z, expectedFraction, 'Each layer uses its own unpooled masked fractions'));
  assert.match(get('minimum-status').textContent, /6 \/ 12/, 'Layer support counts section bins separately instead of pooling measurements');
  app.change('layer-metric', 'n_all'); await app.settle();
  assert.equal(get('minimum').disabled, true);
  app.last().traces.filter(trace => trace.type === 'heatmap').forEach(trace =>
    assert.deepEqual(trace.z, grid(250).n_all.map(row => row.map(Math.log1p)), 'Layer counts stay unfiltered'));

  app.change('layer-view', 'stack'); await app.settle();
  assert.equal(get('camera-controls').hidden, false, 'Exploded stack exposes the same camera controls');
  app.click('camera-top'); await app.settle();
  near(app.last().layout.scene.camera.eye.x, 0, 'Stack supports exact top view');
  assert.ok(app.last().layout.scene.camera.eye.z > 0);

  app.tabs.find(tab => tab.dataset.view === 'channels').onclick(); await app.settle();
  assert.equal(get('camera-controls').hidden, true, '2D channel overlay hides camera controls');
  assert.equal(get('minimum').disabled, true);
  assert.equal(get('bin').disabled, true);
  assert.ok(get('minimum-help').textContent.trim(), 'Fixed-resolution channel view explains disabled cutoff');
  console.log('Explorer controls: camera presets/turns/manual sync, live readout, throttle/flush, denominator masks, count immunity and layered/channel modes passed.');
})().catch(error => { console.error(error); process.exitCode = 1; });
