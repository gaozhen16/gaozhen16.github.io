// TODO: Add doc for class
class Constellation {
    // List of Satellite instances in the constellation
    sats = [];

    // Constructor for Abstract Constellation class
    constructor() {
        // Make sure the class cannot be instantiated
        if (new.target === Constellation) {
            throw new TypeError("Cannot construct \"Constellation\" (Abstract) instances directly");
        }

        // Make sure the class implements these methods
        if (this.updatePositions === undefined) {
            throw new TypeError("Must override method \"updatePositions\"");
        }
    }
}
