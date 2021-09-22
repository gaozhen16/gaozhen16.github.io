// Create Renderer
const renderer = new THREE.WebGLRenderer({antialias:true});
renderer.setSize( window.innerWidth, window.innerHeight );
renderer.setClearColor(0xffffff, 1);
document.body.appendChild( renderer.domElement );

// Create Scene
const scene = new THREE.Scene();

// Create Camera
const camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 100.0 * R_earth );
camera.up.set(0,0,1);
scene.add(camera);

// Orbit controls
const controls = new THREE.OrbitControls( camera, renderer.domElement );
camera.position.y = 4.0 * R_earth;
camera.position.z = 1.0 * R_earth;
camera.lookAt(0,0,0);
controls.update();

// Earth
var _config = {
    radius: R_earth,//半径
    //map: new THREE.TextureLoader().load('https://i.loli.net/2020/11/03/KcXdH9yCf43Y2GZ.jpg'),//加载需要的地球贴图 地球图为宽高 2:1的图 
    map: new THREE.TextureLoader().load('./Files/earth.jpg'),
    //[url=https://sm.ms/image/RXCha2zlgw4j39W][img]https://i.loli.net/2020/11/11/RXCha2zlgw4j39W.jpg[/img][/url]
    //https://i.loli.net/2020/11/11/bwRWgxdQvukV9mf.jpg
}
const earth_geom = new THREE.SphereGeometry(R_earth, 32, 32);
const earth_mat = new THREE.MeshPhongMaterial({
    color: 0xffffff,
    map: _config.map,
});
const earth_mesh = new THREE.Mesh(earth_geom, earth_mat).rotateX(80.1).rotateZ(160.2);
scene.add(earth_mesh);

// Ground Station
const ground_station = new GroundStation(30.0 * Math.PI / 180.0, 0.0);

// Axes of intertial reference frame
const axesHelper = new THREE.AxesHelper(0.5 * R_earth);
scene.add( axesHelper );

// Satellite Constellation
const walker_con = new WalkerConstellation(60.0 * Math.PI / 180.0, 1000, 100, 30, 1.25 * R_earth, scene, ground_station);

// Ambient Light Source
const ambient_light = new THREE.AmbientLight(0xf1f1f1, 1);
scene.add(ambient_light);

// Render loop
var t = 0.0;
function render() {
    t += 10.0;
    requestAnimationFrame(render);
    ground_station.updatePosition(t);
    walker_con.updatePositions(t);

    controls.update();
    renderer.render(scene, camera);
}
render();
