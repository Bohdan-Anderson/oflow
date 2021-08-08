/*jslint sloppy: true, vars: true, plusplus: true, white: true */

var FlowZone = require('./flowZone');

module.exports = FlowCalculator;

/**
 * The heart of the optical flow detection. Implements Lucas-Kande method:
 * http://en.wikipedia.org/wiki/Lucas%E2%80%93Kanade_method
 * Current implementation is not extremely tolerant to garbage collector.
 * This could be improved...
 * 
 * Step dictacts how much of a gradent we use per point
 */
function FlowCalculator(step) {
    this.step = step || 8;
}

FlowCalculator.prototype.calculate = function (oldImage, newImage, width, height) {
    var zones = [];
    var step = this.step;
    var winStep = step * 2 + 1;

    var A2, A1B2, B1, C1, C2;
    var u, v, uu, vv, ww;
    ww = uu = vv = 0;
    var wMax = width - step - 1;
    var hMax = height - step - 1;
    var globalY, globalX, localY, localX;

    var centerX = wMax/2;
    var centreY = hMax/2;

    for (globalY = step + 1; globalY < hMax; globalY += winStep) {
        for (globalX = step + 1; globalX < wMax; globalX += winStep) {
            A2 = A1B2 = B1 = C1 = C2 = 0;
            
            // We create a gradient of the current point
            // Size of the gradient is dictated by step
            for (localY = -step; localY <= step; localY++) {
                for (localX = -step; localX <= step; localX++) {
                    var address = (globalY + localY) * width + globalX + localX;

                    // get the pixel next to +1x -1x +1y -1y
                    var gradX = (newImage[(address - 1) * 4]) - (newImage[(address + 1) * 4]);
                    var gradY = (newImage[(address - width) * 4]) - (newImage[(address + width) * 4]);
                   
                    // Get the previous value for this location (Time)
                    var gradT = (oldImage[address * 4]) - (newImage[address * 4]);

                    A2 += gradX * gradX;
                    A1B2 += gradX * gradY;
                    B1 += gradY * gradY;
                    C2 += gradX * gradT;
                    C1 += gradY * gradT;
                }
            }

            var delta = (A1B2 * A1B2 - A2 * B1);

            if (delta !== 0) {
                /* system is not singular - solving by Kramer method */
                var iDelta = step / delta;
                var deltaX = -(C1 * A1B2 - C2 * B1);
                var deltaY = -(A1B2 * C2 - A2 * C1);

                u = deltaX * iDelta;
                v = deltaY * iDelta;
            } else {
                /* singular system - find optical flow in gradient direction */
                var norm = (A1B2 + A2) * (A1B2 + A2) + (B1 + A1B2) * (B1 + A1B2);
                if (norm !== 0) {
                    var iGradNorm = step / norm;
                    var temp = -(C1 + C2) * iGradNorm;

                    u = (A1B2 + A2) * temp;
                    v = (B1 + A1B2) * temp;
                } else {
                    u = v = 0;
                }
            }

            // Calculate the if it's going towards the middle or going away from the middle
            // for each point - is the line facing towards the center of the frame, or the inverse?
            // once the know the vector what the speed away from the middle is the point moving
            // middle is hMax and wMax
            // added centerX
            // added centreY

            // what is the current distance from centre
            var baseDistance = Math.sqrt(Math.pow(globalX-centerX,2)+Math.pow(globalY-centreY,2))
            var newDistance = Math.sqrt(Math.pow(globalX+u-centerX,2)+Math.pow(globalY+v-centreY,2))
            w = baseDistance-newDistance;


            if (-winStep < u && u < winStep &&
                -winStep < v && v < winStep) {
                uu += u;
                vv += v;
                ww += w;
                zones.push(new FlowZone(globalX, globalY, 0, u, v, w));
            }
        }
    }

    return {
        zones : zones,
        u : uu / zones.length,
        v : vv / zones.length,
        w : ww / zones.length
    };
};
