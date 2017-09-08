// turn on strict mode for supporting browsers
"use strict";

// TODO:
//  - line plots
//  - chords with width (not just chord lines)
//  - tutorials
//  - premade data structures for mouse and human genomes

/**
 * @typedef {Object} CartesianCoords
 * @property {number} x
 * @property {number} y
 */

/**
 * @typedef {Object} PolarCoords
 * @property {number} radians
 * @property {number} distance
 */

/**
 * @typedef {Object} GenoCoords
 * @property {string} chr
 * @property {number} pos
 */

/**
 * @typedef {Object} ValuedGenoCoords
 * @property {string} chr
 * @property {number} pos
 * @property {number} value
 */

/**
 * @typedef {Object} LabeledGenoCoords
 * @property {string} chr
 * @property {number} pos
 * @property {string} text
 */

/**
 * @typedef {Object} Interval
 * @property {number} startPos
 * @property {number} size
 */

/**
 * @typedef {Object} GenoInterval
 * @property {string} chr
 * @property {number} startPos
 * @property {number} size
 */

/**
 * @typedef {Object} ValuedGenoInterval
 * @property {string} chr
 * @property {number} startPos
 * @property {number} size
 * @property {number} value
 */

/**
 * A structure combining three kinds of position coordinates for convenience
 * @typedef {Object} PositionInfo
 * @property {CartesianCoords} cartesian
 * @property {PolarCoords} polar
 * @property {GenoCoords} geno
 */

/**
 * Data structure representing a connection on two points of the genomic coordinate system
 * @typedef {Object} GenoConnection
 * @property {GenoCoords} source
 * @property {GenoCoords} dest
 */

/**
 * determines if a string ends with the given suffix
 * @param {string} str
 * @param {string} suffix
 * @returns {boolean} true iff the given string ends with the given suffix
 */
function endsWith(str, suffix) {
    return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

/**
 * Convert radians to degrees
 */
function radToDeg(rad) {
    return 180.0 * rad / Math.PI;
}

/**
 * Convert degrees to radians
 * @param {number} deg
 * @returns {number} radians
 */
function degToRad(deg) {
    return Math.PI * deg / 180.0;
}

/**
 * Convert polar coordinates to cartesian
 * @param {number} radius
 * @param {number} thetaRad
 * @returns {CartesianCoords}
 */
function polarToCartesian(radius, thetaRad) {
    return {
        "x" : Math.cos(thetaRad) * radius,
        "y" : Math.sin(thetaRad) * radius
    };
}

/**
 * return a value that is the same angle as the given rad but forced to be
 * between the values of 0 and 2pi
 */
function clampRad(rad) {
    rad %= 2.0 * Math.PI;
    if(rad < 0.0) {
        rad += 2.0 * Math.PI;
    }
    return rad;
}

/**
 * Helper function to calculate the exclusive end position (interval.startPos + interval.size)
 * @param {Interval} interval
 * @returns {number} the exclusive end position
 */
function endPosExclu(interval) {
    return interval.startPos + interval.size;
}

/**
 * Helper function to calculate the center position (interval.startPos + interval.size / 2)
 * @param {Interval} interval
 * @returns {number} the center position
 */
function centerPos(interval) {
    return interval.startPos + interval.size / 2;
}

/**
 * this is basically an associated array that makes it easy to track element
 * ordering via the keys member
 * @constructor
 * @param {Array} ordPairs
 *      an array of pairs (that is, an array where each element is a list of
 *      length two). The first element of each pair is a string key and the
 *      second is a value which can be any type.
 */
function OrderedPairs(ordPairs) {
    
    this.length = ordPairs.length;
    
    /** an array containing just the keys in order */
    this.keys = new Array(ordPairs.length);
    
    /** an assoaciative array from the keys to the values */
    this.assocArr = {};
    
    /** an associative array from the keys to their indices */
    this.keyIndices = {};
    
    var i;
    for(i = 0; i < ordPairs.length; i++) {
        var currPair = ordPairs[i];
        var currKey = currPair[0];
        var currVal = currPair[1];

        this.keys[i] = currKey;
        this.keyIndices[currKey] = i;
        this.assocArr[currKey] = currVal;
    }
    
    this.valueAt = function(index) {
        return this.assocArr[this.keys[index]];
    };
}

/**
 * the circular coord system for jaxle
 * @constructor
 * @param {OrderedPairs} chrIntervals
 *      the OrderedPairs object with keys being chromosome names and values
 *      being the {Interval} representing the start position and length
 *      of the chromosome
 * @param {number} spacingRad
 *      how much spacing should we allow in radians between chromosomes
 */
function GenoCoords(chrIntervals, spacingRad) {
    var self = this;

    /**
     * the chrIntervals passed into the constructor
     * @type {OrderedPairs}
     */
    this.chrIntervals = chrIntervals;

    /**
     * an array of d3 scales using the same chromosome ordering as in chrIntervals. The scale domains are in
     * genome coordinates and the range is in radians
     * @type {Array}
     */
    this.genoToRadianScales = new Array(chrIntervals.length);
    var radiansPerUnit = 0;
    
    // create a new scope for init work (avoids polluting this scope
    // with temporary variable names)
    (function() {
        // in order to calc radiansPerUnit we need to accumulate all the chr lengths
        var i;
        var totalAccumLen = 0;
        for(i = 0; i < chrIntervals.length; i++) {
            totalAccumLen += chrIntervals.valueAt(i).size;
        }
        radiansPerUnit =
            (2.0 * Math.PI - spacingRad * chrIntervals.length) /
            totalAccumLen;
        
        // now that we have radians per unit we're ready to calculate all of the
        // linear scales that map from genomic units to radians
        var accumLen = 0;
        for(i = 0; i < chrIntervals.length; i++) {
            
            var currInterval = chrIntervals.valueAt(i);
            var currOffsetRad = spacingRad * i + accumLen * radiansPerUnit;
            
            // rotate 3/4 so that the 1st chromosome starts at the top
            currOffsetRad = clampRad(currOffsetRad + 1.5 * Math.PI);
            
            var currLenRad = radiansPerUnit * currInterval.size;
            
            self.genoToRadianScales[i] =
                d3.scale.linear()
                    .domain([currInterval.startPos, endPosExclu(currInterval)])
                    .range([currOffsetRad, currOffsetRad + currLenRad]);

            accumLen += currInterval.size;
        }
    })();

    /**
     * Converts from the given genome position to radians
     * @param {string} chr the chromosome ID
     * @param {number} pos the position
     * @returns {number} the radians
     */
    this.genoPosToRadians = function(chr, pos) {
        var chrIdx = this.chrIntervals.keyIndices[chr];
        var chrScale = this.genoToRadianScales[chrIdx];
        
        return chrScale(pos);
    };

    /**
     * Converts a radius and genomic position into a cartesian position (where 0,0 is the center of the circle)
     * @param radius the radius
     * @param {string} chr the chromosome ID
     * @param {number} pos the position
     * @returns {CartesianCoords}
     */
    this.genoPosToCartesian = function(radius, chr, pos) {
        var radPos = this.genoPosToRadians(chr, pos);
        return polarToCartesian(radius, radPos);
    };

    /**
     * Converts a position given in radians to a genomic position
     * @param rad the radians to convert to a position
     * @returns {GenoCoords} returns the genomic position or null if there is no valid genome position for the given rad
     */
    this.radiansToGenoPos = function(rad) {
        for(var i = 0; i < this.chrIntervals.length; i++) {
            var chrScale = this.genoToRadianScales[i];
            var startRad = chrScale.range()[0];
            var stopRad = chrScale.range()[1];

            var adj_rad = rad;
            if(adj_rad < startRad) {
                adj_rad += 2.0 * Math.PI;
            }

            if(startRad <= adj_rad && adj_rad <= stopRad) {
                return {
                    'chr': this.chrIntervals.keys[i],
                    'pos': Math.round(chrScale.invert(adj_rad))
                };
            }
        }

        return null;
    };
}

/**
 * A JaxlePlot is the main entry point that allows you to draw circle plots in a given d3 SVG object.
 * @constructor
 * @param {GenoCoords} genoCoords the coordinate system to use for the plot TODO type is wrong
 * @param svg the d3 svg element to the plot to
 */
function JaxlePlot(genoCoords, svg) {
    // keeping track of "this" in javascript hurts my brain so I prefer to work with an explicit "self" here
    var self = this;

    // if rootGroup is defined it means that this object was pre-constructed by self.newGroup() and so we don't have
    // any initialization work to do
    if(typeof(this.rootGroup) === "undefined") {
        (function() {
            var widthStyle = svg.style("width");
            var widthPx = parseInt(widthStyle.substr(widthStyle, widthStyle.length - 2)) / 2;
            var heightStyle = svg.style("height");
            var heightPx = parseInt(heightStyle.substr(heightStyle, heightStyle.length - 2)) / 2;
            self.radius = Math.min(widthPx, heightPx);
            self.rootGroup = svg.append("g").attr("transform", "translate(" + self.radius + ", " + self.radius + ")");
            self.parent = null;
        })();
    }

    /**
     * A convenience function for determining if the given chromosome is part of this plot's genomic coordinate system.
     * @param {string} chr the name of the chromosome
     * @returns {boolean}
     */
    this.chrIsValid = function(chr) {
        return chr in genoCoords.chrIntervals.assocArr;
    };

    /**
     * Creates a new child JaxlePlot that is basically the same as it's parent except that it renders within
     * it's own rootGroup which is a child of the parent's rootGroup. This could be useful for applying styles
     * and transformations to a subset of SVG elements
     * @returns {JaxlePlot} the child plot object
     */
    this.createChildPlot = function() {
        var thisArg = {
            radius: self.radius,
            rootGroup: self.rootGroup.append("g"),
            parent: self
        };
        JaxlePlot.call(thisArg, genoCoords, svg);
        return thisArg;
    };

    // creates an SVG arc path based on the given params
    var _arcPath = function(startRadius, stopRadius, genoInterval) {
        var chrIdx = genoCoords.chrIntervals.keyIndices[genoInterval.chr];
        var chrScale = genoCoords.genoToRadianScales[chrIdx];

        var startRad = chrScale(genoInterval.startPos);
        var stopRad = chrScale(endPosExclu(genoInterval));
        var innerStartPos = polarToCartesian(startRadius, startRad);
        var innerStopPos = polarToCartesian(startRadius, stopRad);
        var outerStartPos = polarToCartesian(stopRadius, stopRad);
        var outerStopPos = polarToCartesian(stopRadius, startRad);

        // Path for an arc:
        // A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        // NOTE: sweep-flag=1 means "positive angle" and 0 means "negative angle"
        return (
            "M " + innerStartPos.x + " " + innerStartPos.y +
            " A " + startRadius + " " + startRadius + " 0 0 1 " + innerStopPos.x + " " + innerStopPos.y +
            " L " + outerStartPos.x + " " + outerStartPos.y +
            " A " + stopRadius + " " + stopRadius + " 0 0 0 " + outerStopPos.x + " " + outerStopPos.y +
            " Z");
    };

    /**
     * Arc post-processing callback
     * @callback arcPostProcFunc
     * @param  svgElem
     *              the d3 object for the SVG path element for the arc
     * @param  {GenoInterval} interval
     *              the arc interval
     * @param  {number} index
     *              index of array element
     */
    /**
     * Draw the given arc shapes
     * @param {number} startRadius
     *          the starting radius for the arc (with respect to the center of the plot)
     * @param {number} radiusExtent
     *          this defines how thick the arc shape will be
     * @param {Array.<GenoInterval>} genoIntervals
     *          all of the genomic intervals to draw arcs for
     * @param {arcPostProcFunc} [postProcFunc]
     *          this callback gives you a chance to modify each arc after it is added to the SVG doc
     * @returns the new d3 group object that all of the arcs are added to. This group will be classed as arc-default
     */
    this.drawArcs = function(startRadius, radiusExtent, genoIntervals, postProcFunc) {
        var fExists = typeof(postProcFunc) !== "undefined";
        var grp = self.rootGroup.append("g").attr("class", "arc-default");

        for(var i = 0; i < genoIntervals.length; i++) {
            if(self.chrIsValid(genoIntervals[i].chr)) {
                var p = grp.append("path").attr("d", _arcPath(startRadius, startRadius + radiusExtent, genoIntervals[i]));
                if(fExists) {
                    postProcFunc(p, genoIntervals[i], i);
                }
            } else {
                console.warn('invalid chromosome in drawArcs: ' + genoIntervals[i].chr);
            }
        }

        return grp;
    };

    var _arcLinePath = function(radiusPx, genoInterval) {
        var chrIdx = genoCoords.chrIntervals.keyIndices[genoInterval.chr];
        var chrScale = genoCoords.genoToRadianScales[chrIdx];

        var startRad = chrScale(genoInterval.startPos);
        var stopRad = chrScale(endPosExclu(genoInterval));
        var startPos = polarToCartesian(radiusPx, startRad);
        var stopPos = polarToCartesian(radiusPx, stopRad);

        // Path for an arc:
        // A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        // NOTE: sweep-flag=1 means "positive angle" and 0 means "negative angle"
        return (
            "M " + startPos.x + " " + startPos.y +
            " A " + radiusPx + " " + radiusPx + " 0 0 1 " + stopPos.x + " " + stopPos.y);
    };

    /**
     * Arc post-processing callback
     * @callback labelPostProcFunc
     * @param  svgElem
     *              the d3 object for the SVG text element
     * @param  {LabeledGenoCoords} labeledCoords
     *              the text coords
     * @param  {number} index
     *              index of array element
     */
    /**
     * Draw labels on the plot
     * @param {number} params.radiusPx the radius
     * @param {Array.<LabeledGenoCoords>} params.labels all of the labels (NOTE: these objects will be modified)
     * @param {labelPostProcFunc} [params.postProcFunc] post processing function for label elements
     * @return the d3 group containing all of the labels. This group is classed as genomic-label-default
     */
    this.drawGenomicLabels = function(params) {
        var fExists = typeof(params.postProcFunc) !== 'undefined';
        var grp = self.rootGroup.append('g').attr('class', 'genomic-label-default');

//        // We will need to create and merge label clusters in order to avoid overlapping labels. This
//        // means that we have to order labels on a radian scale from greatest to least on a radian scale,
//        // however we don't want to take the naive radian comparison approach because a cluster can span
//        // the zero boundary.
//        //
//        //
//        var labelClusters = [];
//
//        function mergeLabelClusters(clust1, clust2) {
//            var i1 = 0;
//            var i2 = 0;
//
//            while(i1 < clust1.labeledCoords.length || i2 < clust2.labeledCoords.length) {
//                if(i1 >= clust1.labeledCoords.length)
//            }
//        }
//
        params.labels.forEach(function(label, i) {
            if(self.chrIsValid(label.chr)) {
//                //var p = grp.append("path").attr("d", _arcLinePath(radiusPx, genoIntervals[i]));
//                var labelPos = genoCoords.genoPosToCartesian(params.radiusPx, label.chr, label.pos);
//                var radians = genoCoords.genoPosToRadians(label.chr, label.pos);
//                var txt = grp.append('text')
//                    .text(label.text)
//                    .style('dominant-baseline', 'central');
//                var txtBBox = txt.node().getBBox();
//                var txtHeight = txtBBox.height;
//                //var txtWidth = txtBBox.width; // TODO may need to use this for scaling
//                var sizeRadians = Math.tan(txtHeight / params.radiusPx);
//
//                var currCluster = {
//                    sizeRadians: sizeRadians,
//                    centerRadians: radians,
//                    startRadians: clampRad(radians - sizeRadians / 2),
//                    stopRadians: clampRad(radians + sizeRadians / 2)//,
//                    //OK Keith
//                };
//
                var txtTrans = null;
                var deg = radToDeg(radians);
                if(deg > 90 && deg < 270){
                    txtTrans =
                        "translate(" + labelPos.x + "," + labelPos.y +
                        ")rotate(" + (deg + 180) + ")";
                    txt.style('text-anchor', 'end');
                } else {
                    txtTrans =
                        "translate(" + labelPos.x + "," + labelPos.y +
                        ")rotate(" + deg + ")";
                }
                txt.attr("transform", txtTrans);

                if(fExists) {
                    params.postProcFunc(txt, label, i);
                }
            } else {
                console.warn('invalid chromosome in drawLabels: ' + label.chr);
            }
        });

        return grp;
    };

    /**
     * Draw arcs (part of the circumference) as a line.
     * @param {number} radiusPx the radius of the circle to draw an arc for
     * @param {Array.<GenoInterval>} genoIntervals
     * @param {arcPostProcFunc} [postProcFunc]
     *          this callback gives you a chance to modify each arc line after it is added to the SVG doc
     * @return the d3 SVG group containing all of the arc lines. This group is classed as arc-line-default.
     */
    this.drawArcLines = function(radiusPx, genoIntervals, postProcFunc) {
        var fExists = typeof(postProcFunc) !== "undefined";
        var grp = self.rootGroup.append("g").attr("class", "arc-line-default");

        for(var i = 0; i < genoIntervals.length; i++) {
            if(self.chrIsValid(genoIntervals[i].chr)) {
                var p = grp.append("path").attr("d", _arcLinePath(radiusPx, genoIntervals[i]));
                if(fExists) {
                    postProcFunc(p, genoIntervals[i], i);
                }
            } else {
                console.warn('invalid chromosome in drawArcLines: ' + genoIntervals[i].chr);
            }
        }

        return grp;
    };

    var _chordPath = function(raduisPx, startChr, startPos, stopChr, stopPos) {
        var startPosPx = genoCoords.genoPosToCartesian(raduisPx, startChr, startPos);
        var stopPosPx = genoCoords.genoPosToCartesian(raduisPx, stopChr, stopPos);

        return (
            "M " + startPosPx.x + "," + startPosPx.y +
            " Q 0,0 " + stopPosPx.x + "," + stopPosPx.y);
    };

    /**
     * Chord line post-processing callback
     * @callback chordPostProcFunc
     * @param  svgElem
     *              the d3 object for the SVG path element for the chord
     * @param  {GenoConnection} data
     *              the chord data
     * @param  {number} index
     *              index of array element
     */
    /**
     * Draw chord lines (chord lines connect two locations on the plot crossing through the middle of the plot)
     * @param {number} radiusPx
     *          the radius
     * @param {Array.<GenoConnection>} chordLines
     *          the chord lines to draw
     * @param {chordPostProcFunc} [postProcFunc]
     * @return the d3 SVG group containing all of the chord lines. This group is classed as chord-line-default.
     */
    this.drawChordLines = function(radiusPx, chordLines, postProcFunc) {
        var fExists = typeof(postProcFunc) !== "undefined";
        var grp = self.rootGroup.append("g").attr("class", "chord-line-default");

        for(var i = 0; i < chordLines.length; i++) {

            var chordLine = chordLines[i];
            var path = _chordPath(
                radiusPx,
                chordLine.source.chr,
                chordLine.source.pos,
                chordLine.dest.chr,
                chordLine.dest.pos
            );
            var p = grp.append("path").attr("d", path);
            if(fExists) {
                postProcFunc(p, chordLine, i);
            }
        }

        return grp;
    };

    /**
     * Take a chromosome and position and return a formatted string value
     * @callback chrPosFormatCallback
     * @param {string} chr the chromosome ID
     * @param {number} pos the position
     */
    /**
     * Draw ticks as specified by params
     * @param {number} params.radius
     *          the radius where the tick marks start
     * @param {number} params.tickInterval
     *          the interval between ticks (in genomic distance). If there are minor
     *          ticks this is the interval between minor ticks.
     * @param {number} params.tickSize
     *          the tick size (in pixels). If there are minor ticks this is the major tick size.
     * @param {number} [params.ticksPerLabel]
     *          the number of ticks rendered for every label (default is to not render labels)
     * @param {chrPosFormatCallback} [params.labelFormatFunc]
     *          the default behavior is to return the pos value converted to a string
     * @param {number} [params.labelTickMargin=2]
     *          the distance between the end of a major tick and the start of a label
     * @param {number} [params.ticksPerMajorTick=1]
     *          the total number of ticks to draw for every major tick (eg. a value
     *          of indicates that every tick is a major tick while a value of 2 indicates that every other tick is
     *          a major tick)
     * @param {number} [params.minorTickSize=params.tickSize/2.0]
     *          the minor tick size (in pixels)
     * @param [params.chrs=genoCoords.chrIntervals.keys]
     *          defines which chromosomes ticks should be drawn for
     * @return the new d3 group object that all of the ticks are added to.
     *         This group will be classed as tick-default. There are also subgroups
     *         for tick labels classed as tick-label-default, minor ticks classed as minor-tick-default
     *         and major ticks classed as major-tick-default
     */
    this.drawTicks = function(params) {
        if(typeof(params.ticksPerMajorTick) === "undefined") {
            params.ticksPerMajorTick = 1;
        }
        if(typeof(params.minorTickSize) === "undefined") {
            params.minorTickSize = params.tickSize / 2.0;
        }
        if(typeof(params.labelFormatFunc) === "undefined") {
            params.labelFormatFunc = function(chr, pos) {return pos.toString()};
        }
        if(typeof(params.labelTickMargin) === "undefined") {
            params.labelTickMargin = 2;
        }
        if(typeof(params.chrs) === "undefined") {
            params.chrs = genoCoords.chrIntervals.keys;
        }
        var renderLabels = typeof(params.ticksPerLabel) !== "undefined";

        // create groups in such a way that we can organize the SVG classes
        var grp = self.rootGroup.append("g");
        var tickGroup = grp.append("g").attr("class", "tick-default");
        var majorTickGroup = tickGroup.append("g").attr("class", "major-tick-default");
        var minorTickGroup;
        if(params.ticksPerMajorTick > 1) {
            minorTickGroup = tickGroup.append("g").attr("class", "minor-tick-default");
        }
        var tickLabelGroup;
        if(renderLabels) {
            tickLabelGroup = grp.append("g").attr("class", "tick-label-default");
        }

        params.chrs.forEach(function(chr) {
            var chrInterval = genoCoords.chrIntervals.assocArr[chr];
            var currPos = chrInterval.startPos;

            // adjust currPos so that ticks are offset from 0
            if(currPos % params.tickInterval !== 0) {
                currPos += params.tickInterval - (currPos % params.tickInterval);
            }

            // iterate through all of the tick positions
            var tickCount = 0;
            while(currPos < endPosExclu(chrInterval)) {
                // set size and group according to if this is a major or minor tick
                var currSize;
                var currGrp;
                if(tickCount % params.ticksPerMajorTick === 0) {
                    // it's a major tick
                    currGrp = majorTickGroup;
                    currSize = params.tickSize;
                } else {
                    // it's a minor tick
                    currGrp = minorTickGroup;
                    currSize = params.minorTickSize;
                }

                // draw the tick
                var tickStart = genoCoords.genoPosToCartesian(params.radius, chr, currPos);
                var tickStop = genoCoords.genoPosToCartesian(params.radius + currSize, chr, currPos);
                currGrp.append("line")
                    .attr("x1", tickStart.x)
                    .attr("y1", tickStart.y)
                    .attr("x2", tickStop.x)
                    .attr("y2", tickStop.y);

                // draw label
                if(renderLabels && tickCount % params.ticksPerLabel === 0) {
                    var labelRadius = params.radius + params.tickSize + params.labelTickMargin;
                    var labelPos = genoCoords.genoPosToCartesian(labelRadius, chr, currPos);

                    var radians = genoCoords.genoPosToRadians(chr, currPos);
                    var txtTrans =
                        "translate(" + labelPos.x + "," + labelPos.y +
                        ")rotate(" + (90 + radToDeg(radians)) + ")";
                    tickLabelGroup.append("text")
                        .text(params.labelFormatFunc(chr, currPos))
                        .attr("text-anchor", "middle")
                        .attr("transform", txtTrans);
                }

                // increment position and count
                currPos += params.tickInterval;
                tickCount++;
            }
        });

        return grp;
    };

    /**
     * function that takes a chromosome string and returns the label string to use. The default behavior
     * is to use the chromosome ID directly
     * @callback chrFormatCallback
     * @param {string} chr the chromosome ID
     * @returns {string}
     */
    /**
     * Draw chromosome labels
     * @param {number} params.radius
     *          the radius where the labels will be rendered
     * @param [params.chrs=genoCoords.chrIntervals.keys]
     *          defines which chromosomes labels should be rendered for
     * @param {chrFormatCallback} [params.labelFormatFunc]
     *          function that takes a chromosome string and returns the label string to use. The default behavior
     *          is to use the chromosome ID directly
     * @return the new d3 group object that all of the labels are added to.
     *         This group will be classed as chr-label-default
     */
    this.drawChrLabels = function(params) {
        if(typeof(params.labelTickMargin) === "undefined") {
            params.labelTickMargin = 2;
        }
        if(typeof(params.chrs) === "undefined") {
            params.chrs = genoCoords.chrIntervals.keys;
        }
        if(typeof(params.labelFormatFunc) === "undefined") {
            params.labelFormatFunc = function(chr){return chr;};
        }

        var grp = self.rootGroup.append("g").attr("class", "chr-label-default");

        params.chrs.forEach(function(chr) {
            var chrInterval = genoCoords.chrIntervals.assocArr[chr];
            var chrMidPoint = (chrInterval.startPos + chrInterval.size) / 2.0;
            var labelPos = genoCoords.genoPosToCartesian(params.radius, chr, chrMidPoint);
            var radians = genoCoords.genoPosToRadians(chr, chrMidPoint);
            var txtTrans =
                "translate(" + labelPos.x + "," + labelPos.y +
                ")rotate(" + (90 + radToDeg(radians)) + ")";
            grp.append("text")
                .text(params.labelFormatFunc(chr))
                .attr("text-anchor", "middle")
                .attr("transform", txtTrans);
        });

        return grp;
    };

    /**
     * Histogram bar post-processing callback
     * @callback histPostProcFunc
     * @param  svgElem
     *              the d3 object for the SVG path element for a histogram bar
     * @param  {ValuedGenoInterval} data
     *              the histogram bar data
     * @param  {number} index
     *              index of array element
     */
    /**
     * Draw a histogram with the following data
     * @param params.scale
     *          a d3 scale object where the domain is given in the units provided by the values in the
     *          params.data array and the range is in radius units
     * @param {Array.<ValuedGenoInterval>} params.data
     *          an array of intervals and values each of which will be rendered as a bar on the histogram
     * @param {histPostProcFunc} [params.postProcFunc]
     *          this callback gives you a chance to modify each histogram bar after it is added to the SVG doc
     * @return the new d3 group object that all of the histogram bars are added to.
     *         This group will be classed as hist-bar-default
     */
    this.drawHistogram = function(params) {
        var scale = params.scale;
        var hasPostProcFunc = typeof(params.postProcFunc) !== "undefined";
        var startRadius = typeof(params.baseline) === "undefined" ? scale.range()[0] : scale(params.baseline);

        var grp = self.rootGroup.append("g").attr("class", "hist-bar-default");

        for(var i = 0; i < params.data.length; i++) {
            var valuedInterval = params.data[i];

            if(self.chrIsValid(valuedInterval.chr)) {
                var stopRadius = scale(valuedInterval.value);
                var p = grp.append("path").attr("d", _arcPath(startRadius, stopRadius, valuedInterval));
                if(hasPostProcFunc) {
                    params.postProcFunc(p, valuedInterval, i);
                }
            } else {
                console.warn('invalid chromosome in drawHistogram: ' + valuedInterval.chr);
            }
        }

        return grp;
    };

    this.drawLineGraph = function(params) {
        var scale = params.scale;
        var data = params.data;
        var hasPostProcFunc = typeof(params.postProcFunc) !== "undefined";

        var grp = self.rootGroup.append("g").attr("class", "line-graph-default");

        if(data.length >= 2) {
            var lineFunction = d3.svg.line()
                .x(function(d) { return d.x; })
                .y(function(d) { return d.y; })
                .interpolate("linear");

            var lastBreakIndex = 0;
            var currIndex = 1;
            var renderContigSegments = function() {
                var dataSubset = data.slice(lastBreakIndex, currIndex);
                var dataSubsetXY = dataSubset.map(function(elem) {
                    var radiusPx = scale(elem.value);
                    return genoCoords.genoPosToCartesian(radiusPx, elem.chr, elem.pos);
                });
                grp.append("path").attr("d", lineFunction(dataSubsetXY));

                lastBreakIndex = currIndex;
            };

            for(; currIndex < data.length; currIndex++) {
                var valuedInterval = data[currIndex];
                var chrChanged = valuedInterval.chr !== data[currIndex - 1].chr;
                if(chrChanged) {
                    renderContigSegments();
                }
            }
            renderContigSegments();
        }

        return grp;
    };

    /**
     * Uses the d3.mouse function to get the current mouse event and uses the mouse event's xy position to calculate
     * cartesian, polar and genomic coordinates
     * @returns {PositionInfo}
     */
    this.mousePositionInfo = function() {
        var mouseXY = d3.mouse(svg.node());
        var x = mouseXY[0] - self.radius;
        var y = mouseXY[1] - self.radius;
        var radians = clampRad(Math.atan2(y, x));

        return {
            cartesian: {x: x, y: y},
            polar: {
                radians: clampRad(Math.atan2(y, x)),
                distance: Math.sqrt(x * x + y * y)
            },
            geno: genoCoords.radiansToGenoPos(radians)
        };
    };
}
