import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { OcctKernel } from 'occt-wasm';
import { registerKernel, OcctWasmAdapter, mesh, toBufferGeometryData } from 'brepjs';
import { BevelGear } from 'brep/bevelGear';

// --- brep-js setup
const kernel = await OcctKernel.init();
registerKernel('occt-wasm', OcctWasmAdapter.fromKernel(kernel));
// ---

// A meshing pair of spiral bevel gears on perpendicular axes.
const mod = 2;
const z1 = 8;
const z2 = 21;
const gamma1 = Math.atan2(z1, z2);
const gamma2 = Math.PI / 2 - gamma1;
const beta = Math.PI / 6;

const pinion = new BevelGear({
  numberOfTeeth: z1,
  module: mod,
  height: 10,
  coneAngle: 2 * gamma1,
  helixAngle: beta,
});

const wheel = new BevelGear({
  numberOfTeeth: z2,
  module: mod,
  height: 10,
  coneAngle: 2 * gamma2,
  helixAngle: -beta,
});

// place the pinion against the wheel (meshing position is baked into the solid)
pinion.meshTo(wheel);

const pinionSolid = pinion.buildSolid();
const wheelSolid = wheel.buildSolid();

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x16181d);

const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.1, 1000);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

function createThreeMesh(solid: any) {
  const meshData = mesh(solid);
  const geoData = toBufferGeometryData(meshData);
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(geoData.position, 3));
  geometry.setAttribute('normal', new THREE.BufferAttribute(geoData.normal, 3));
  if (geoData.index) {
    geometry.setIndex(new THREE.BufferAttribute(geoData.index, 1));
  }
  const material = new THREE.MeshStandardMaterial({
    color: 0xb0b0b0,
    metalness: 0.25,
    roughness: 0.55,
  });
  return new THREE.Mesh(geometry, material);
}

scene.add(createThreeMesh(wheelSolid));
scene.add(createThreeMesh(pinionSolid));

scene.add(new THREE.HemisphereLight(0xffffff, 0x33363f, 1.2));
const keyLight = new THREE.DirectionalLight(0xffffff, 1.6);
keyLight.position.set(40, -60, 50);
scene.add(keyLight);
const fillLight = new THREE.DirectionalLight(0xffffff, 0.5);
fillLight.position.set(-50, 40, -20);
scene.add(fillLight);

// Frame the pair.
const bounds = new THREE.Box3().setFromObject(scene);
const sphere = bounds.getBoundingSphere(new THREE.Sphere());
const center = sphere.center.clone();
center.z -= sphere.radius * 0.12; // most of the mass sits low, around the wheel
camera.up.set(0, 0, 1);
camera.position.copy(center).add(
  new THREE.Vector3(0.8, -1, 0.45).normalize().multiplyScalar(sphere.radius * 2),
);
camera.near = sphere.radius / 50;
camera.far = sphere.radius * 20;
camera.updateProjectionMatrix();

const controls = new OrbitControls(camera, renderer.domElement);
controls.target.copy(center);
controls.enableDamping = true;

window.addEventListener('resize', () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});

renderer.setAnimationLoop(() => {
  controls.update();
  renderer.render(scene, camera);
});
