(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.oflow = f()}})(function(){var define,module,exports;return (function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
/*global window,  */

var FlowCalculator = require('./flowCalculator.js');

module.exports = CanvasFlow;

/**
 * A high level interface to capture optical flow from the <canvas> tag.
 * The API is symmetrical to webcamFlow.js
 *
 * Usage example:
 *  var flow = new VideoFlow();
 *
 *  // Every time when optical flow is calculated
 *  // call the passed in callback:
 *  flow.onCalculated(function (direction) {
 *      // direction is an object which describes current flow:
 *      // direction.u, direction.v {floats} general flow vector
 *      // direction.zones {Array} is a collection of flowZones.
 *      //  Each flow zone describes optical flow direction inside of it.
 *  });
 *  // Starts capturing the flow from webcamera:
 *  flow.startCapture();
 *  // once you are done capturing call
 *  flow.stopCapture();
 */
function CanvasFlow(defaultCanvasTag, zoneSize) {
    var calculatedCallbacks = [],
        canvas = defaultCanvasTag,
        ctx,
        width,
        height,
        oldImage,
        loopId,
        calculator = new FlowCalculator(zoneSize || 8),

        requestAnimFrame = window.requestAnimationFrame       ||
                           window.webkitRequestAnimationFrame ||
                           window.mozRequestAnimationFrame    ||
                           window.oRequestAnimationFrame      ||
                           window.msRequestAnimationFrame     ||
                           function( callback ) { window.setTimeout(callback, 1000 / 60); },
        cancelAnimFrame =  window.cancelAnimationFrame ||
                           window.mozCancelAnimationFrame,
        isCapturing = false,

        getCurrentPixels = function () {
            return ctx.getImageData(0, 0, width, height).data;
        },
        calculate = function () {
            var newImage = getCurrentPixels();
            if (oldImage && newImage) {
                var zones = calculator.calculate(oldImage, newImage, width, height);
                calculatedCallbacks.forEach(function (callback) {
                    callback(zones);
                });
            }
            oldImage = newImage;
        },

        initView = function () {
            width = canvas.width;
            height = canvas.height;
            ctx = canvas.getContext('2d');
        },
        animloop = function () {
            if (isCapturing) {
                loopId = requestAnimFrame(animloop);
                calculate();
            }
        };

    if (!defaultCanvasTag) {
        var err = new Error();
        err.message = "Video tag is required";
        throw err;
    }

    this.startCapture = function () {
        // todo: error?
        isCapturing = true;
        initView();
        animloop();
    };
    this.stopCapture = function () {
        cancelAnimFrame(loopId);
        isCapturing = false;
    };
    this.onCalculated = function (callback) {
        calculatedCallbacks.push(callback);
    };
    this.getWidth = function () { return width; };
    this.getHeight = function () { return height; };
}

},{"./flowCalculator.js":2}],2:[function(require,module,exports){
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

},{"./flowZone":3}],3:[function(require,module,exports){
module.exports = FlowZone;

function FlowZone(x, y, z, u, v, w) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.u = u;
    this.v = v;
    this.w = w;
}

},{}],4:[function(require,module,exports){
module.exports = {
  WebCamFlow: require('./webcamFlow'),
  VideoFlow: require('./videoFlow'),
  CanvasFlow: require('./canvasFlow'),
  FlowZone: require('./flowZone'),
  FlowCalculator: require('./flowCalculator')
};

},{"./canvasFlow":1,"./flowCalculator":2,"./flowZone":3,"./videoFlow":5,"./webcamFlow":6}],5:[function(require,module,exports){
/*global window */

var FlowCalculator = require('./flowCalculator');
module.exports = VideoFlow;

/**
 * A high level interface to capture optical flow from the <video> tag.
 * The API is symmetrical to webcamFlow.js
 *
 * Usage example:
 *  var flow = new VideoFlow();
 *
 *  // Every time when optical flow is calculated
 *  // call the passed in callback:
 *  flow.onCalculated(function (direction) {
 *      // direction is an object which describes current flow:
 *      // direction.u, direction.v {floats} general flow vector
 *      // direction.zones {Array} is a collection of flowZones.
 *      //  Each flow zone describes optical flow direction inside of it.
 *  });
 *  // Starts capturing the flow from webcamera:
 *  flow.startCapture();
 *  // once you are done capturing call
 *  flow.stopCapture();
 */
function VideoFlow(defaultVideoTag, zoneSize) {
    var calculatedCallbacks = [],
        canvas,
        video = defaultVideoTag,
        ctx,
        width,
        height,
        oldImage,
        loopId,
        calculator = new FlowCalculator(zoneSize || 8),

        requestAnimFrame = window.requestAnimationFrame       ||
                           window.webkitRequestAnimationFrame ||
                           window.mozRequestAnimationFrame    ||
                           window.oRequestAnimationFrame      ||
                           window.msRequestAnimationFrame     ||
                           function( callback ) { window.setTimeout(callback, 1000 / 60); },
        cancelAnimFrame =  window.cancelAnimationFrame ||
                           window.mozCancelAnimationFrame,
        isCapturing = false,

        getCurrentPixels = function () {
            width = video.videoWidth;
            height = video.videoHeight;
            canvas.width  = width;
            canvas.height = height;

            if (width && height) {
                ctx.drawImage(video, 0, 0);
                var imgd = ctx.getImageData(0, 0, width, height);
                return imgd.data;
            }
        },
        calculate = function () {
            var newImage = getCurrentPixels();
            if (oldImage && newImage) {
                var zones = calculator.calculate(oldImage, newImage, width, height);
                calculatedCallbacks.forEach(function (callback) {
                    callback(zones);
                });
            }
            oldImage = newImage;
        },

        initView = function () {
            width = video.videoWidth;
            height = video.videoHeight;

            if (!canvas) { canvas = window.document.createElement('canvas'); }
            ctx = canvas.getContext('2d');
        },
        animloop = function () {
            if (isCapturing) {
                loopId = requestAnimFrame(animloop);
                calculate();
            }
        };

    if (!defaultVideoTag) {
        var err = new Error();
        err.message = "Video tag is required";
        throw err;
    }

    this.startCapture = function () {
        // todo: error?
        isCapturing = true;
        initView();
        animloop();
    };
    this.stopCapture = function () {
        cancelAnimFrame(loopId);
        isCapturing = false;
    };
    this.onCalculated = function (callback) {
        calculatedCallbacks.push(callback);
    };
    this.getWidth = function () { return width; };
    this.getHeight = function () { return height; };
}

},{"./flowCalculator":2}],6:[function(require,module,exports){
/*global navigator, window */

var VideoFlow = require('./videoFlow');
module.exports = WebCamFlow;

/**
 * A high level interface to capture optical flow from the web camera.
 * @param defaultVideoTag {DOMElement} optional reference to <video> tag
 *   where web camera output should be rendered. If parameter is not
 *   present a new invisible <video> tag is created.
 * @param zoneSize {int} optional size of a flow zone in pixels. 8 by default
 * @param cameraFacing {string} optional direction camera is facing (either
 * 'user' or 'environment') used to give preference to a particular mobile
 * camera. If matching camera is not found, any available one will be used.
 *
 * Usage example:
 *  var flow = new WebCamFlow();
 *
 *  // Every time when optical flow is calculated
 *  // call the passed in callback:
 *  flow.onCalculated(function (direction) {
 *      // direction is an object which describes current flow:
 *      // direction.u, direction.v {floats} general flow vector
 *      // direction.zones {Array} is a collection of flowZones.
 *      //  Each flow zone describes optical flow direction inside of it.
 *  });
 *  // Starts capturing the flow from webcamera:
 *  flow.startCapture();
 *  // once you are done capturing call
 *  flow.stopCapture();
 */
function WebCamFlow(defaultVideoTag, zoneSize, cameraFacing, onFail) {
    var videoTag,
        isCapturing,
        localStream,
        calculatedCallbacks = [],
        selectedVideoSource,
        desiredDevice,
        videoFlow,
        onWebCamFail = function(e) {
            if(e.name === "NotAllowedError"){
                window.alert('You have denied access to your camera.');
            } else {
                window.alert('getUserMedia() is not supported in your browser.');
            }
            if (onFail) {
                onFail();
            }
        },
        onWebCamSucceed = function(stream) {
            isCapturing = true;
            localStream = stream;
            if ("srcObject" in videoTag) {
                videoTag.srcObject = stream;
            } else {
                videoTag.src = window.URL.createObjectURL(stream);
            }
            if (stream) {
                videoTag.play();
                videoFlow.startCapture(videoTag);
                videoFlow.onCalculated(gotFlow);
            }
        },
        gotFlow = function(direction) {
            calculatedCallbacks.forEach(function (callback) {
                callback(direction);
            });
        },
        initCapture = function() {
            if (!videoFlow) {
                videoTag = defaultVideoTag || window.document.createElement('video');
                videoTag.setAttribute('autoplay', true);
                videoFlow = new VideoFlow(videoTag, zoneSize);
            }

            if (window.MediaStreamTrack.getSources) {
                window.MediaStreamTrack.getSources(function(sourceInfos) {
                    for (var i = 0; i < sourceInfos.length; i++) {
                        if (sourceInfos[i].kind === 'video' && confirm(sourceInfos[i])){
                            selectedVideoSource = sourceInfos[i].id;
                            // if camera facing requested direction is found, stop search
                            if (sourceInfos[i].facing === cameraFacing) {
                                break;
                            }
                            break;
                        }
                    }

                    desiredDevice = { optional: [{sourceId: selectedVideoSource}] };
                    navigator.mediaDevices.getUserMedia({ video: desiredDevice })
                        .then(onWebCamSucceed)
                        .catch(onWebCamFail);
                });
            } else if(navigator.mediaDevices.enumerateDevices) {
                navigator.mediaDevices.enumerateDevices().then(
                    function(sourceInfos){
                        for (var i = 0; i < sourceInfos.length; i++) {
                            if(sourceInfos[i].kind == "videoinput" && confirm(sourceInfos[i].label)){
                                selectedVideoSource = sourceInfos[i].deviceId;
                                break;
                            }
                        }
                        
                        desiredDevice = { optional: [{sourceId: selectedVideoSource}] };
                        navigator.mediaDevices.getUserMedia({ video: desiredDevice })
                            .then(onWebCamSucceed)
                            .catch(onWebCamFail);
                    }
                );
            }
        }

    // our public API
    this.startCapture = function () {
        if (!isCapturing) {
            initCapture();
        }
    };
    this.onCalculated = function (callback) {
        calculatedCallbacks.push(callback);
    };
    this.stopCapture = function() {
        isCapturing = false;
        if (videoFlow) { videoFlow.stopCapture(); }
        if (videoTag) { videoTag.pause(); }
        if (localStream) { localStream.stop(); }
    };
}

},{"./videoFlow":5}]},{},[4])(4)
});
