// TODO: Add doc for this class
class WalkerConstellation extends Constellation {
    // Constructor for the WalkerConstellation class
    constructor(i, T, P, F, a, scene, ground_station) {
        // Call Constellation constructor first
        super();

        // Initialize the planes of satellites
        for (var plane = 0; plane < P; plane++) {
            // Compute the RAAN of the satellites in the plane
            const Omega = plane * 2.0 * Math.PI / P;

            // Initialize the satellites in the plane
            for (var n = 0; n < T / P; n++) {
                // Compute the initial argument of periapsis of the satellite
                const nu_0 = (n * 2.0 * Math.PI * P / T) + (plane * 2.0 * Math.PI * F / T);

                // Initialize the satellite
                this.sats.push(new Satellite(a, 0.0, i, 0.0, Omega, nu_0, scene, ground_station));
            }
        }
    }

    // TODO: Add doc for function
    updatePositions(t) {
        for (var i = 0; i < this.sats.length; i++) {
            this.sats[i].updatePosition(t);
        }
    }
}
