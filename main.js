var CLIP_ID = 'clip-states';
var graphMapsOn = 1;
var boundingBoxCoords;
var w = 1280,
    h = 800;
var globalScale =700;

var projection = d3.geo.azimuthal()
    .mode("equidistant")
    .origin([-98, 38])
    .scale(700)
    .translate([640, 360]);

var path = d3.geo.path()
    .projection(projection);

var svg = d3.select("body").insert("svg:svg", "h2")
    .attr("width", w)
    .attr("height", h);

var states = svg.append("svg:g")
    .attr("id", "states");

var cells = svg.append("svg:g")
    .attr("id", "cells");

var svg_ex = d3.select("#svgcontainer");

var renderedRegions = svg.append("svg:g")
    .attr("id", "rendered-regions");

var line = d3.svg.line();
var lineCardinal = d3.svg.line()
    .x(function(d) { return d.X; })
    .y(function(d) { return d.Y; }).
    interpolate("cardinal-closed");


function polygon(ring) {
  var polygon = [ring];
  ring.push(ring[0]); // add closing coordinate
  //if (d3.geo.area({type: "Polygon", coordinates: polygon}) > 2 * Math.PI) ring.reverse(); // fix winding order
  return polygon;
}


function constructDefaultPoints(defaultLevels,stateLatLon,week){

  var reportByState = {};
  Object.keys(defaultLevels).forEach(function(regionName){

    var region = defaultLevels[regionName];
    region.states.forEach(function(state){

      var defaultLatLon = stateLatLon.filter(function(defaultStateLatLon){
        return defaultStateLatLon.state === state;
      }).pop();

      var aReport = reportByState[state] = {
        lat :  parseFloat(defaultLatLon.latitude),
        lng : parseFloat(defaultLatLon.longitude),
        report : {}
      };

      if(!region.weeks) return;

      Object.keys(region.weeks).forEach(function(subject){
        var dateRange = region.weeks[subject];
        var level;
        if(week < dateRange.onset){
          level = 0;
        }else if(week >= dateRange.onset && week < dateRange.peakStart){
          level = 1;
        }else if(week >= dateRange.peakStart && week < dateRange.peakEnd){
          level = 3;
        }else if(week >= dateRange.peakEnd && week < dateRange.end){
          level = 1;
        }else{
          level = 0;
        }

        aReport.report[subject] = {level : level};
      });
    });
  });

  return reportByState;
}


d3.json("resources/default-levels.json", function(defaultLevels) {
  d3.json("resources/us.json", function(us) {

      var selectedStates = {type: "GeometryCollection", geometries: us.objects.states.geometries},
          selectionBoundary = topojson.mesh(us, selectedStates, function(a, b) { return a === b; }),
          selection = {type: "MultiPolygon", coordinates: selectionBoundary.coordinates.map(polygon)};
          bounds = d3.geom.polygon(selectionBoundary);

      svg.append("svg:defs")
          .append('svg:clipPath')
          .attr('id',CLIP_ID)
          .append('path')
          .datum(selection)
          //.attr("class", "state selected")      //can re-enable for debugging
          .attr("d", path);

    d3.csv("resources/state_latlon.csv", function(stateLatLon) {

      var week = 0; //moment().week();

      window.setInterval(function(){
        week = (week + 1) % 52;
	
	var div = document.getElementById('week-title');

	div.innerHTML = 'Week ' + week;

        var reportByState = constructDefaultPoints(defaultLevels,stateLatLon,week);

        render(reportByState,stateLatLon);
      }, 1000);
    });
  });
});



function render(reportByState,stateLatLon){

    function fill(d,i){
      var state = stateLatLon[i].state;
      var doc = reportByState[state];
      if(!doc){ 
        return 0;
      }

      var level = doc.report[compound].level;
      switch(level){
        case 0 :
          return 'green';
        case 1 : 
          return 'yellow';
        case 2 : 
          return 'orange';
        case 3 : 
          return 'red';
      }
    }

    var positions = [];

    stateLatLon.forEach(function(stateAvgLatLon) {
      positions.push(projection([+stateAvgLatLon.longitude, +stateAvgLatLon.latitude]));
    });

    // Compute the Voronoi diagram of state' projected positions.
    var polygons = d3.geom.voronoi(positions);//.map(function(cell) { return selectionBoundary.clip(cell); });
    var compound = 'grasses';

    //-----johns method calls------
    generateBoundingBox(globalScale);
    findConnections(polygons,stateLatLon, reportByState, compound);
    convertIndexToPolygons(polygons);
    //    removeDuplicatePoints();
    //    createHull();
    regions = transformRegionsToXYFormat(regions);
    polygonsCopy = transformPolygonsToXYFormat(polygonsCopy);
    roundToZeroPolygons();
    mergePolygons();
    boundingBox();
    if (graphMapsOn) graphRegions();
    cleanUp();
    

    var g = cells.selectAll("g")
        .data(stateLatLon);

    var enter = g.enter()

    enter.append("svg:g")
        .append("svg:path")
	.attr("d", function(d, i) { return line(polygons[i]); })
        .attr("class", "cell")
        .attr('clip-path','url(#' + CLIP_ID + ')');
    enter.append("svg:circle")
        .attr("cx", function(d, i) { return positions[i][0]; })
        .attr("cy", function(d, i) { return positions[i][1]; })
        .attr("r", 1.5);
    cells.selectAll("g")
        .data(stateLatLon)
        .attr('fill',fill);

}

//regions starts empty, then gets filled with connected polygons, which then get converted to points, then duplicate points removed,  then a copy of regions is made called regionsHull, then a hull is calculated for all the regions and stored in regionsHull.  This data is graphed.
var regions = [
	       [		], //level 0
	       [		], //level 1
	       [		], //level 2
	       [		] //level 3
	       ],
    distanceArray = [
		     [],
		     [],
		     [],
		     []
		     ],
    connectednessGraph = [],
    mergeStatus = [],
    globalRenderMergeCount = 0,
    polygonsCopy;


function cleanUp () { //reset the regions and regionsHull arrays
    polygonsCopy = null;
    mergeStatus = [];
    connectednessGraph = [];
    regions = [
	       [		], //level 0
	       [		], //level 1
	       [		], //level 2
	       [		] //level 3
	       ];
    distanceArray = [
	       [		], //level 0
	       [		], //level 1
	       [		], //level 2
	       [		] //level 3
		     ];
    globalRenderMergeCount =0;
}

function generateBoundingBox(scale) {
    var normal = 700; //the bounding box is normalized according to a scale of 700
    boundingBoxCoords = [[
    {X: (350*scale)/normal, Y: (150*scale)/normal },
    {X: (950*scale)/normal, Y:(150*scale)/normal },
    {X: (950*scale)/normal , Y:(550*scale)/normal},
    {X: (350*scale)/normal, Y:(550*scale)/normal }
			  ]];


}

function createHull () { //takes regions array and computes hull in regionsHull for each connected region
    regionsHull  = $.extend(true, [], regions); //deep copy of regions array
    var firstPoint;
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    regionsHull[level][connectedRegion] = d3.geom.hull(regionsHull[level][connectedRegion]);
	}
    }
}

function subtractPoly (subj_paths, clip_paths) {
    if (subj_paths ==undefined || clip_paths==undefined) return;
    subj_paths = [subj_paths]; //add external array so that it works with jsclipper
    clip_paths = [clip_paths]; // //add external array so that it works with jsclipper
    var scale = 1;
    ClipperLib.JS.ScaleUpPaths(subj_paths, scale);
    ClipperLib.JS.ScaleUpPaths(clip_paths, scale);
    var cpr = new ClipperLib.Clipper();
    cpr.AddPaths(subj_paths, ClipperLib.PolyType.ptSubject, true);
    cpr.AddPaths(clip_paths, ClipperLib.PolyType.ptClip, true);
    var subject_fillType = ClipperLib.PolyFillType.pftNonZero;
    var clip_fillType = ClipperLib.PolyFillType.pftNonZero;
    var clipTypes = [ClipperLib.ClipType.ctDifference];
    var solution_paths = new ClipperLib.Paths();
    cpr.Execute(clipTypes[0], solution_paths, subject_fillType, clip_fillType);
    return solution_paths[0];
}

function addPoly (subj_paths, clip_paths) {
    if (subj_paths ==undefined || clip_paths==undefined) return;
    var scale = 1;
    ClipperLib.JS.ScaleUpPaths(subj_paths, scale);
    ClipperLib.JS.ScaleUpPaths(clip_paths, scale);
    var cpr = new ClipperLib.Clipper();
    cpr.AddPaths(subj_paths, ClipperLib.PolyType.ptSubject, true);
    cpr.AddPaths(clip_paths, ClipperLib.PolyType.ptClip, true);
    var subject_fillType = ClipperLib.PolyFillType.pftNonZero;
    var clip_fillType = ClipperLib.PolyFillType.pftNonZero;
    var clipTypes = [ClipperLib.ClipType.ctUnion];
    var solution_paths = new ClipperLib.Paths();
    cpr.Execute(clipTypes[0], solution_paths, subject_fillType, clip_fillType);
    return solution_paths;
}


function removeClosePoints() { //removes points that are too close
    var limit = 15;  //anything closer than "limit" pixels gets removed
    //go through the arrays, removing the 2nd point that corresponds to a line distance < limit
    for (var level = 0; level < distanceArray.length; level++) {
	if (distanceArray[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < distanceArray[level].length; connectedRegion++) {
	    if (distanceArray[level][connectedRegion].length==0) continue;
	    for (var regionArea =0 ; regionArea < regions[level][connectedRegion].length; regionArea++) {
		if (distanceArray[level][connectedRegion][regionArea].length==0) continue;
		for (var distanceIndex = distanceArray[level][connectedRegion][regionArea].length -1 ; distanceIndex >= 0; distanceIndex--) { //we work backwards since we'll be deleting
		    if (distanceArray[level][connectedRegion][regionArea][distanceIndex]<limit) {
			if (distanceIndex!=(distanceArray[level][connectedRegion][regionArea].length-1)) {
			    //delete the 2nd point
			    if (regions[level][connectedRegion][regionArea][distanceIndex+1])
				regions[level][connectedRegion][regionArea].splice(distanceIndex+1,1); 
			    //we remove it if it exists...might not if we moved the last point beforehand
			}
			else {
			    //we're at the last index in hte distance array, which corresopnds to distnace between last and first point in connectedregions array
			    regions[level][connectedRegion][regionArea].splice(distanceIndex,1);
			    //we remove THE FIRST POINT (not second), which corresponds to last point in connectedregions array
			}
		    }
		}
	    }
	}
    }
}

function lineDistance( point1, point2) {
    var xs =0;
    var ys = 0;

    xs = point2.X - point1.X;
    xs = xs* xs;

    ys = point2.Y - point1.Y;
    ys = ys* ys;

    return Math.sqrt(xs + ys);
}

function distanceBetweenPoints() { //[0] corresponds to distance between [0] and [1] in connectedRegionPoints 
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion]) distanceArray[level].push([]);//push an empty arry to correspond to connectedRegion
	    if (regions[level][connectedRegion].length==0) continue;
	    for (var regionArea =0 ; regionArea < regions[level][connectedRegion].length; regionArea++) {
		if (regions[level][connectedRegion][regionArea]) distanceArray[level][connectedRegion].push([]);//push an empty arry to correspond to connectedRegion
		for (var regionPoint = 0; regionPoint < regions[level][connectedRegion][regionArea].length; regionPoint++) { //we start at 1
		    if (regionPoint!=(regions[level][connectedRegion][regionArea].length-1)) {
			var dist = lineDistance(regions[level][connectedRegion][regionArea][regionPoint], regions[level][connectedRegion][regionArea][regionPoint+1]);
			distanceArray[level][connectedRegion][regionArea].push(dist);
		    }
		    else {
			//case for last and first point distance
			var dist = lineDistance(regions[level][connectedRegion][regionArea][regionPoint], regions[level][connectedRegion][regionArea][0]);
			distanceArray[level][connectedRegion][regionArea].push(dist);
		    }
		}
	    }
	}
    }
}


function polygonAreaHelper(X, Y, numPoints) 
{ 
    area = 0;         // Accumulates area in the loop
    j = numPoints-1;  // The last vertex is the 'previous' one to the first

    for (i=0; i<numPoints; i++)
	{ area = area +  (X[j]+X[i]) * (Y[j]-Y[i]); 
	    j = i;  //j is previous vertex to i
	}
    return area/2;
}

function polygonAreaMain(polygon) {
    if (polygon[0] instanceof Array) polygon=polygon[0]; 
    //this takes care of the fact that jsclip will return inner array
    var xPts=[],yPts=[];
    for (var i = 0; i < polygon.length; i++) {
	xPts.push(polygon[i].X);
	yPts.push(polygon[i].Y);
    }
    return Math.abs(polygonAreaHelper(xPts,yPts, yPts.length));
}

function isGrowing(numBefore, numAfter) {
    return (numAfter >=numBefore) || numBefore==undefined;
}



function generateMergeStatus(polygons) {
    for (var i=0; i<polygons.length; i++) {
	connectednessGraph.push([]);
    }
}

function roundToZero() { //rounds polygon points to zero decimal point so it's compatible with jsclipper
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    for (var regionPoly = 0; regionPoly < regions[level][connectedRegion].length; regionPoly++) { //we start at 1n
		polygonA = regions[level][connectedRegion][regionPoly][0];
		for (var y = 0; y< polygonA.length; y++) {
		    //this takes care of the fact that jsclip will return inner array
		    var aX = polygonA[y].X.toFixed(0),//rounded to zero decimals b/c that's what jsclipperdoes
			aY = polygonA[y].Y.toFixed(0);
		    polygonA[y].X = +aX;
		    polygonA[y].Y = +aY;
		    

		}
		regions[level][connectedRegion][regionPoly]=polygonA;
	    }
	}
    }
}

function roundToZeroPolygons() { 
    //rounds polygon points to zero decimal point so it's compatible with jsclipper
    for (var i = 0; i < polygonsCopy.length; i++) {
	//	for (var j = 0 ; j <polygonsAry[i].length; j++){
	polygonA = polygonsCopy[i];
	for (var y = 0; y< polygonA.length; y++) {
	    var aX = polygonA[y].X.toFixed(0),//rounded to zero decimals b/c that's what jsclipperdoes
		aY = polygonA[y].Y.toFixed(0);
	    polygonA[y].X = +aX;
	    polygonA[y].Y = +aY;
	}
    }
}

function mergePolygons() { //for given connected region, merges polygon points
    mergeStatus.length=connectednessGraph.length;  
    roundToZero();
    roundToZeroPolygons();
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    //	    var firstPolyIndex = regionsIndex[level][connectedRegion][0];
	    var firstPoly = [regions[level][connectedRegion][0]];//***IMPORTANT add array around to match jsclip
	    var firstPolyIndex = regionsIndex[level][connectedRegion][0];
	    mergeStatus[firstPolyIndex]=true;

	    growingPoly = merge(firstPoly, regionsIndex[level][connectedRegion][1]);
	    regions[level][connectedRegion]=growingPoly;
	    //	    renderAfterMerge($.extend(true, [], growingPoly)); //deep copy of regions arraygrowingPolygon);
	}
    }

        distanceBetweenPoints();
        removeClosePoints();
}

function boundingBox(){ //goes through regions and cuts around the US MAP
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    var subj_paths = regions[level][connectedRegion];
	    var clip_paths = boundingBoxCoords;
	    var scale = 1;
	    ClipperLib.JS.ScaleUpPaths(subj_paths, scale);
	    ClipperLib.JS.ScaleUpPaths(clip_paths, scale);
	    var cpr = new ClipperLib.Clipper();
	    cpr.AddPaths(subj_paths, ClipperLib.PolyType.ptSubject, true);
	    cpr.AddPaths(clip_paths, ClipperLib.PolyType.ptClip, true);
	    var subject_fillType = ClipperLib.PolyFillType.pftNonZero;
	    var clip_fillType = ClipperLib.PolyFillType.pftNonZero;
	    var clipTypes = [ClipperLib.ClipType.ctIntersection];
	    var solution_paths = new ClipperLib.Paths();
	    cpr.Execute(clipTypes[0], solution_paths, subject_fillType, clip_fillType);
	    regions[level][connectedRegion]=solution_paths;
	}
    }
}

function merge(growingPoly, currentPolyIndex) {
    if (mergeStatus[currentPolyIndex])  {
	return growingPoly;}
    else {
        //MERGE
	if (!polygonsCopy[currentPolyIndex]) return growingPoly; //this somehow fixed an edge case i don't understand
        growingPoly = addPoly(growingPoly, [polygonsCopy[currentPolyIndex]]); //***IMPORTNAT adding external array for jsclipper
	//renderAfterMerge($.extend(true, [], growingPoly)); //deep copy of regions arraygrowingPolygon);
        mergeStatus[currentPolyIndex] = true;
        if (connectednessGraph[currentPolyIndex].length > 0) {
            for (var a = 0; a < connectednessGraph[currentPolyIndex].length; a++) {
		var polyIndex = connectednessGraph[a];
                if (mergeStatus[polyIndex]) continue;
		growingPoly = merge(growingPoly, connectednessGraph[currentPolyIndex][a]);
            }
	}
	return growingPoly;
    }
}

function renderAfterMerge(growingPolygon) { //test render after each merge to know what hte fuck is going on.
    //    var svg_ex, cont = document.getElementById('svgcontainer');
    var pathString = paths2string(growingPolygon, 1);
    svg_ex.append("svg:svg")
	.attr("width", 1600)
	.attr("height", 1600)
	.append("svg:g")
	.append("path").data(growingPolygon)
	.attr("id", "rendered-regions-" + globalRenderMergeCount)
	.attr("fill","blue")
	.attr("d" , /*pathString);*/ function (d, i) { 
	  return lineCardinal(growingPolygon[0]) +"Z "/* +
							       lineCardinal(growingPolygon[1]) + "Z" }*/});
    globalRenderMergeCount++;
}

function isDone(allocArray) {
    for (var a = 0; a < allocArray.length; a++) {
	if (!allocArray[a]) return false;
    }
    return true;
}

function removeDuplicatePoints () { // for a given connected region, removes points that appear twice
   for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    for (var regionPoint = 0; regionPoint < regions[level][connectedRegion].length; regionPoint++) {
		if (regions[level][connectedRegion][regionPoint].length==0) continue; // this should never happen, but just in case
		var currentPoint = regions[level][connectedRegion][regionPoint];
		for (var a =0; a < regions[level][connectedRegion].length; a++) {
		    var comparePoint = regions[level][connectedRegion][a];
		    if (currentPoint[0] === comparePoint[0] && currentPoint[1] === comparePoint[1]) {
			//we have duplicate. delete
			regions[level][connectedRegion].splice(a,1); 
			a--; // this cancels out a++ and means i stay in place after deleting in line prev
		    }
		}
	    }
	}
    }
}

var regionsIndex;

//this takes the regions array of indices, and replaces it with copies of the actual polygons
function convertIndexToPolygons(polygons) { 
    regionsIndex  = $.extend(true, [], regions); //deep copy of regions array, which at the time contains indices
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regions[level].length; connectedRegion++) {
	    if (regions[level][connectedRegion].length==0) continue;
	    for (var regionPoly = 0; regionPoly < regions[level][connectedRegion].length; regionPoly++) {
		if (regions[level][connectedRegion][regionPoly].length==0) continue; // this should never happen, but just in case
		var indexOfPoly = regions[level][connectedRegion][regionPoly];
		var polygonCopy = polygons[indexOfPoly];
		regions[level][connectedRegion][regionPoly] = polygonCopy;
	    }
	}
    }
}


function generateCnxnGraph(polygons) {
    for (var i=0; i<polygons.length; i++) {
	connectednessGraph.push([]);
    }
}

function alreadyConnected(polygonIndexA, polygonIndexB) {
    for (var a = 0; a< connectednessGraph[polygonIndexA].length; a++ ) {
	if (connectednessGraph[polygonIndexA][a]===polygonIndexB) return true;
    }
    return false;
}


function findConnections (polygons,stateLatLon, reportByState, compound) {
    generateCnxnGraph(polygons);
    polygonsCopy  = $.extend(true, [], polygons); //deep copy of polygons, for use in merge
    for (var i  = 1; i < polygons.length; i++) {
	var state = stateLatLon[i].state;
	var doc = reportByState[state];
	if(!doc) continue;
	var level = doc.report[compound].level;
	var thisRegionIndex = isRegion(i,level); // returns >=0 if in region
	for (var j = 1; j < polygons.length; j++) {
	    if (i==j) continue;  //same polygon
	    var compareState = stateLatLon[j].state;
	    var compareDoc = reportByState[compareState];
	    if(!compareDoc) continue;
	    var compareLevel = compareDoc.report[compound].level;
	    if (level === compareLevel && isConnected(polygons[i],polygons[j])) {
		if (!alreadyConnected(i,j)) {
		    connectednessGraph[i].push(j);
		    connectednessGraph[j].push(i); }//THATS ALL!
		var existingRegionIndex = isRegion(j,level); //check if match is in existing region
		if (existingRegionIndex >=0  &&  thisRegionIndex==-1)  { //case 1: he is in region, and i am not, so add me
		    regions[level][existingRegionIndex].push(i);
		    thisRegionIndex = existingRegionIndex; 
		}
		else if (existingRegionIndex>=0  &&  thisRegionIndex >=0 ) {//case 2:  he is in region, and so am i. merge us.
		    if (existingRegionIndex == thisRegionIndex) continue; //they're in the same region
		    regions[level][existingRegionIndex] = 
			regions[level][existingRegionIndex].concat(regions[level][thisRegionIndex]);
		    regions[level].splice(thisRegionIndex, 1); //delete this region
		    if (thisRegionIndex < existingRegionIndex) //this means that when I deleted at this RegionIndex, existingRegionIndex ->-1
			thisRegionIndex = existingRegionIndex - 1;
		    else
			thisRegionIndex = existingRegionIndex;
		}
		else if (existingRegionIndex == -1 && thisRegionIndex>=0) { //case 3: i'm in a region but he isn't.  add him
		    regions[level][thisRegionIndex].push(j);
		} 
		else { //case 4: neither of us are in a region, create one and add both of us
		    thisRegionIndex = regions[level].push([i,j]) - 1; // .push returns new length of array,
		}
	    }
	    if (thisRegionIndex===-1 && j===(polygons.length - 1)) {
		regions[level].push([i]);
	    }
	}
    }
}


function isConnected (polygonA, polygonB) {
    for (var y = 0; y < polygonA.length; y++) {  //compare points
	for (var z=0; z< polygonB.length; z++) {
	    //IMPORTANT.  rounding to two decimal places b/c of errors with matching points
	    var aX = polygonA[y][0].toFixed(2),
		bX = polygonB[z][0].toFixed(2),
		aY = polygonA[y][1].toFixed(2),
		bY = polygonB[z][1].toFixed(2); 

	    // reassigning rounded point values in polygons for future use
	    polygonA[y][0] = +aX;
	    polygonB[z][0] = +bX;
	    polygonA[y][1] = +aY;
	    polygonB[z][1] = +bY;
	    
	    if (aX===bX && aY===bY) //points match
		return true;
	}
    }
    return false;
}

//this is an alternative written for JSCLipper format
function isConnectedXY (polygonA, polygonB) {
    if (polygonA[0] instanceof Array) polygonA=polygonA[0]; //this takes care of the fact that jsclip will return inner array
    if (polygonB[0] instanceof Array) polygonB=polygonB[0];
    for (var y = 0; y < polygonA.length; y++) {  //compare points
	for (var z=0; z< polygonB.length; z++) {
	    //IMPORTANT.  rounding to two decimal places b/c of errors with matching points
	    var aX = polygonA[y].X.toFixed(0),//rounded to zero decimals b/c that's what jsclipperdoes
		bX = polygonB[z].X.toFixed(0),
		aY = polygonA[y].Y.toFixed(0),
		bY = polygonB[z].Y.toFixed(0);
	    polygonA[y].X = +aX;
	    polygonB[z].X = +bX;
	    polygonA[y].Y = +aY;
	    polygonB[z].Y = +bY;
	    // reassigning rounded point values in polygons for future use
	    if (aX===bX && aY===bY) //points match
		return true;
	}
    }
    return false;
}

/* Checks of a polygon (given by its index) has already been allocated to a region */
function isRegion(index, level) {
    for (var x = 0 ; x < regions[level].length; x++) {
	for (var y=0; y < regions[level][x].length; y++) {
	    if (regions[level][x][y]==index) return x; //regions[level][x]; //return the region array index
	}
    }
    return -1; //failed to find it in an existing region
}

function graphRegions() { //each connected region is drawn
    for (var level = 0; level < regions.length; level++) {
	if (regions[level].length==0) continue;
	for (var m = 0; m < regions[level].length; m++) {
	    if (regions[level][m].length==0) continue;
	    renderedRegions.append("svg:g")
		.append("path").data(regions[level][m])
		.attr("fill-rule","evenodd")
		.attr("d", function(d) { 
			return lineCardinal(d) + "Z";
		    })
		.attr("class", "level-"+level)
		.attr('clip-path','url(#' + CLIP_ID + ')');}
    }
}



function transformRegionsToXYFormat (regionsAry) {
    for (var level = 0; level < regionsAry.length; level++) {
	if (regionsAry[level].length==0) continue;
	for (var connectedRegion = 0; connectedRegion < regionsAry[level].length; connectedRegion++) {
	    if (regionsAry[level][connectedRegion].length==0) continue;
	    for (var regionPoly =0; regionPoly< regionsAry[level][connectedRegion].length; regionPoly++) {
		regionsAry[level][connectedRegion][regionPoly] = arrayToObjArray(regionsAry[level][connectedRegion][regionPoly]);
		regionsAry[level][connectedRegion][regionPoly] = [regionsAry[level][connectedRegion][regionPoly]];
		// i added an array around this so that everything's in jsclipper format
	    }
	}
    }
    return regionsAry;
}



function transformPolygonsToXYFormat (polygonsAry) {
    for (var i = 0; i < polygonsAry.length; i++) {
	polygonsAry[i] = arrayToObjArray(polygonsAry[i]);
    }
    return polygonsAry;
}

function arrayToObjArray(a) {
    var b = [];

    for (var x = 0; x < a.length; x++) {
	b.push({"X": a[x][0], "Y":a[x][1]});
    }
    return b;
}
