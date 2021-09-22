// Earth rotation rate
const earth_rot_rate = 2.0 * Math.PI / 86164.0905; // rad / s

// TODO: Add doc for class
class GroundStation {
    // Constructor for the GroundStation
    constructor(lat, lon) {
        // Save the latitude and longitude
        this.lat = lat;
        this.lon = lon;

        // GroundStation geometry in renderiing
        const geometry = new THREE.SphereGeometry(0.001 * R_earth, 32, 32);

        // GroundSation material in rendering
        const material = new THREE.MeshLambertMaterial( { color: 0x000000 } );

        // GroundSation mesh in rendering
        this.mesh = new THREE.Mesh(geometry, material);

        // Add the GroundStation mesh to the scene
        scene.add(this.mesh);
    }

    // TODO: Add doc for this function
    position(t) {
        // Compute longitude shifted by earth rotation in earth inertial frame
        const lambda = (this.lon + earth_rot_rate * t) % (2.0 * Math.PI);

        // Compute the position in the earth centered inertial frame
        const r_ = new THREE.Vector3(Math.cos(lambda) * Math.cos(this.lat),
                                     Math.sin(lambda) * Math.cos(this.lat),
                                     Math.sin(this.lat));
        r_.multiplyScalar(R_earth);

        // Return the position vector
        return r_;
    }

    // TODO: Add doc for this function
    updatePosition(t) {
        // Update the position of the ground station
        this.mesh.position.copy(this.position(t));
    }
}
