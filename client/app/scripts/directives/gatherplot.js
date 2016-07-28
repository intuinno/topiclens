(function() {
    'use strict';

    angular.module('myApp.directives')
        .directive('gatherplot', ['$http', '$log',
                function($http, $log) {
                    return {
                        restrict: 'EA',
                        scope: {
                            data: "=",
                            config: "=",
                            border: "=",
                            round: "=",
                            xdim: "@",
                            ydim: "@",
                            shapeRenderingMode: "=",
                            dimsum: "=",
                            context: "=",
                            comment: "=",
                            onClick: '&'
                        },

                        link: function(scope, iElement, iAttrs) {

                            //Constants and Setting Environment variables 

                            var margin = 80;
                            var zoom;


                            var maxDotSize = 4;

                            if (scope.config.matrixMode === true) {
                                margin = 5;
                                maxDotSize = 5;
                            }

                            var socket = io.connect('http://0.0.0.0:5004/subtopic');
                            var opt = { epsilon: 20 };
                            var tsne = new subtsnejs.subtSNE(opt);
                            var tsne_animation;
                            var mouseOverTimer;
                            var Y;
                            var subCluster;
                            var subText;

                            var shiftDowned = false;

                            var xRange, yRange;



                            var width = 1040;
                            var height = 820;
                            var outerWidth = width + 2 * margin;
                            var outerHeight = height + 2 * margin;
                            var colorNominal = d3.scale.category20();
                            var color;
                            var colorScaleForHeatMap = d3.scale.linear()
                                .range(["#98c8fd", "08306b"])
                                .interpolate(d3.interpolateHsl);
                            var renderData;

                            var xValue, yValue; //Function for getting value of X,Y position 
                            var xOriginalValue, yOriginalValue;
                            var xScale, yScale;
                            var xAxis, yAxis;
                            var xMap, yMap;
                            var nest = {};

                            var defaultBinSize = 10;

                            var marginForBorderOfAxis = 0.5; //Margin for Border Of Axis


                            var marginClusterRatio = 0.1; //Ratio of margin in the cluster 
                            var node;

                            var dimSetting = {};

                            var svg, svgGroup, nodeGroup, brushGroup, xAxisNodes, yAxisNodes;
                            var maskGroup;
                            var tooltip;
                            var clusterControlBox;

                            var labelDiv;

                            scope.config.binSiz = defaultBinSize;

                            var initialLensSize = 100;

                            var initialRectLensWidth = 100;

                            var initialRectLensHeight = 100;

                            var initialCircleLensR = 50;

                            var initialHistLensWidth = 120;
                            var initialHistLensHeight = 90;

                            var initialInnerRadiusOfPieLens = 20;

                            var brush = d3.svg.brush();
                            var shiftKey;


                            var lensInfo = {};
                            // dimsum = {};


                            var initializeSVG = function() {


                                d3.select("body")
                                    .attr("tabindex", 1)
                                    .on("keydown.brush", keyflip)
                                    .on("keyup.brush", keyflip)
                                    .each(function() {
                                        this.focus();
                                    });

                                // .value("title");

                                if (!scope.config.matrixMode) {

                                    labelDiv = d3.select(iElement[0])
                                        .append("div")
                                        .attr("class", "btn-group")
                                        .html('<a class="btn btn-default" title="Pan and Zoom" id="toolbarPanZoom"><i class="fa fa-search-plus"></i></a><a class="btn btn-default" title="Select" id="toolbarSelect"><i class="fa fa-square-o"></i></a><a class="btn btn-default" title="Reset" id="toolbarReset"><i class="fa fa-undo"></i></a>');
                                }
                                svg = d3.select(iElement[0])
                                    .append("svg:svg");

                                svgGroup = svg.append("g")
                                    .attr("transform", "translate(" + margin + "," + margin + ")");

                                maskGroup = svgGroup.append("g")
                                    .attr("class", "masks");

                                nodeGroup = maskGroup.append("g")
                                    .attr("class", "nodes");

                                nodeGroup.append("rect")
                                    .attr("class", "overlay")
                                    .attr("width", width)
                                    .attr("height", height);

                                brushGroup = svg.append("g")
                                    .attr("transform", "translate(" + margin + "," + margin + ")");


                                xAxisNodes = svgGroup.append("g")
                                    .attr("class", "x axis")
                                    .attr("transform", "translate(0," + height + ")");


                                yAxisNodes = svgGroup.append("g")
                                    .attr("class", "y axis");

                                tooltip = d3.select("body").append("div")
                                    .attr("class", "tooltip")
                                    .style("opacity", 0)
                                    .style('width','auto')
                                    .style('height','auto')
                                    .style('border','1px solid black');

                                clusterControlBox = d3.select("body").append("div")
                                    .attr("class", "clusterControlBox")
                                    .style("opacity", 0);



                            };

                            initializeSVG();

                            // on window resize, re-render d3 canvas
                            window.onresize = function() {
                                return scope.$apply();
                            };

                            scope.$watch(function() {
                                return angular.element(window)[0].innerWidth;
                            }, function() {
                                return scope.handleConfigChange(renderData, scope.config);
                            });

                            // watch for data changes and re-render
                            scope.$watch('data', function(newVals, oldVals) {
                                return scope.renderDataChange(newVals, scope.config);

                            }, false);

                            // watch for Config changes and re-render

                            scope.$watch('config', function(newVals, oldVals) {
                                // debugger;
                                return scope.handleConfigChange(renderData, newVals);
                            }, true);

                            scope.$watch(function() {
                                return scope.border;
                            }, function(newVals, oldVals) {
                                return scope.renderBorderChange(newVals);
                            }, false);

                            scope.$watch(function() {
                                return scope.round;
                            }, function(newVals, oldVals) {
                                return scope.renderRoundChange(newVals);
                            }, false);

                            scope.$watch(function() {
                                return scope.comment;
                            }, function(newVals, oldVals) {
                                if (newVals === true) {
                                    return scope.handleConfigChange(renderData, scope.config);
                                }
                            }, false);

                            scope.$watch(function() {
                                return scope.shapeRenderingMode;
                            }, function(newVals, oldVals) {
                                return scope.renderShapeRenderingChange(newVals);
                            }, true);


                            scope.$watch(function() {
                                return scope.dimsum;
                            }, function(newVals, oldVals) {
                                return scope.handleDimsumChange(newVals);

                            }, true);



                            scope.handleDimsumChange = function(newDimsum) {

                                if (!node) {
                                    return;
                                }

                                if (!scope.dimsum) {

                                    return;
                                }

                                if (!scope.dimsum.selectionSpace) {

                                    return;
                                }



                                node.classed("selected", function(d) {

                                    if (scope.dimsum.selectionSpace.indexOf(d.id) == -1) {

                                        d.selected = false;
                                    } else {

                                        d.selected = true;
                                    }

                                    return d.selected;
                                });

                                scope.dimsum.source = "gatherplot";

                            };


                            scope.renderBorderChange = function(isBorder) {

                                svgGroup.selectAll(".dot")
                                    .style("stroke", function(d) {
                                        return isBorder ? 'black' : 'none';
                                    });

                            };

                            scope.renderRoundChange = function(isRound) {

                                svgGroup.selectAll(".dot")
                                    .transition()
                                    .duration(500)
                                    .attr("rx", function(d) {
                                        return isRound ? +d.nodeWidth / 2 : 0;
                                    })
                                    .attr("ry", function(d) {
                                        return isRound ? +d.nodeWidth / 2 : 0;
                                    });

                            };

                            scope.renderShapeRenderingChange = function(newShapeRendering) {

                                svgGroup.selectAll(".dot")
                                    .style("shape-rendering", newShapeRendering);

                            };

                            var reloadDataToSVG = function() {

                                svgGroup.selectAll("*").remove();

                                maskGroup = svgGroup.append("g")
                                    .attr("class", "masks");

                                nodeGroup = maskGroup.append("g")
                                    .attr("class", "nodes");

                                maskGroup.selectAll("rect").remove();

                                drawBackground();

                                if (scope.config.matrixMode === false) {

                                    node = nodeGroup.selectAll(".dot")
                                        .data(scope.data)
                                        .enter().append("rect")
                                        .attr("class", "dot")
                                        .attr('fill-opacity',0.8)
                                        .on("mouseover", function(d) {
                                            var pageX = d3.event.pageX,
                                                pageY = d3.event.pageY;
                                            mouseOverTimer = setTimeout(function() {
                                                nodeGroup.selectAll('text').style('opacity',0);

                                                var authorList = d.authors.split(', ');
                                                var tooltipString = '<b>Authors:<br/>';
                                                for(var i=0;i<authorList.length;i++) {
                                                    if(i==4) tooltipString+='<br/>'+authorList[i];
                                                    else tooltipString+=authorList[i];
                                                }
                                                tooltipString+='<br/>Title:<br/>';
                                                var titleWords = d.title.split(' ');
                                                var accumulatedLength=0;
                                                var maxLineLength = 40;
                                                for(var i=0;i<titleWords.length;i++) {
                                                    accumulatedLength+=titleWords[0].length+1;
                                                    if(accumulatedLength>maxLineLength) {
                                                        accumulatedLength=0;
                                                        tooltipString+='<br/>'+titleWords[i]+' ';
                                                    }
                                                    else tooltipString+=titleWords[i]+' ';
                                                }
                                                tooltipString+='</b>'

                                                tooltip.transition()
                                                    .duration(100)
                                                    .style("opacity", 0.9);


                                                //tooltip.html(scope.xdim + ":" + xOriginalValue(d) + "<br/>" + scope.ydim + ":" + yOriginalValue(d) + "<br/>" + scope.config.colorDim + ":" + colorOriginalValue(d) + "")
                                                tooltip.html(tooltipString)
                                                    .style("left", (pageX + 5) + "px")
                                                    .style("top", (pageY - 28) + "px")
                                                    .style('background-color','white');

                                            }, 700);

                                        })
                                        .on("mouseout", function(d) {
                                            clearTimeout(mouseOverTimer);
                                            nodeGroup.selectAll('text').style('opacity',1);
                                            tooltip.style("opacity", 0);
                                        })
                                        .on("mousedown", function(d) {
                                            if (d3.event.shiftKey) d3.select(this).classed("selected", d.selected = !d.selected);
                                            else node.classed("selected", function(p) {
                                                return p.selected = d === p;
                                            });
                                        });

                                } else {

                                    nodeGroup.selectAll(".dot")
                                        .data(scope.data)
                                        .enter().append("rect")
                                        .attr("class", "dot");

                                    svg.on("mouseover", function(d) {
                                            tooltip.transition()
                                                .duration(200)
                                                .style("opacity", 0.9);


                                            tooltip.html("<h3>" + scope.xdim + " vs " + scope.ydim + "</h3>")
                                                .style("left", (d3.event.pageX + 5) + "px")
                                                .style("top", (d3.event.pageY - 28) + "px");
                                        })
                                        .on("mouseout", function(d) {
                                            tooltip.transition()
                                                .duration(500)
                                                .style("opacity", 0);
                                        })
                                        .on("click", function(d) {

                                            return scope.onClick({
                                                item: {
                                                    xDim: scope.xdim,
                                                    yDim: scope.ydim
                                                }
                                            });
                                        });


                                }


                                dimSetting = {};


                            };

                            var identifyAndUpdateDimDataType = function() {

                                for (var i = 0; i < scope.config.dims.length; i++) {

                                    var dim = scope.config.dims[i];
                                    dimSetting[dim] = {};
                                    dimSetting[dim].dimType = identifyDimDataType(dim);
                                    prepareDimSettingKeys(dim);

                                }

                            };

                            var prepareDimSettingKeys = function(dim) {

                                var currentDimSetting = dimSetting[dim];

                                if (currentDimSetting.dimType === 'ordinal') {

                                    //doBinningAndSetKeys(dim);
                                    currentDimSetting.isBinned = true;


                                } else {

                                    setKeysFromOriginalData(dim);
                                    currentDimSetting.isBinned = false;

                                }


                            };


                            var doBinningAndSetKeys = function(dimName, numBin) {

                                var currentDimSetting = dimSetting[dimName];

                                currentDimSetting.binnedData = scope.data.map(binningFunc(dimName, numBin));

                            };

                            var binningFunc = function(dimName, numBin) {

                                var minValue = d3.min(scope.data, function(d) {
                                    return +d[dimName];
                                });
                                var maxValue = d3.max(scope.data, function(d) {
                                    return +d[dimName];
                                });

                                var encodingBinScale = d3.scale.linear()
                                    .range([0, numBin - 1])
                                    .domain([minValue, maxValue]);

                                var decodingBinScale = d3.scale.linear()
                                    .domain([0, numBin - 1])
                                    .range([minValue, maxValue]);

                                var binKeys = d3.range(0, numBin, 1);

                                binKeys = binKeys.map(function(d) {
                                    return decodingBinScale(d + 0.5);
                                });


                                dimSetting[dimName].halfOfBinDistance = (decodingBinScale(1) - decodingBinScale(0)) / 2;

                                dimSetting[dimName].keyValue = initializeKeyValueObject(binKeys);

                                return function(d) {

                                    return decodingBinScale(Math.floor(encodingBinScale(d[dimName])) + 0.5);
                                };

                            };

                            var setKeysFromOriginalData = function(dim) {

                                if (!dim) {

                                    return '';

                                }

                                var nest = d3.nest()
                                    .key(function(d) {
                                        return d[dim];
                                    })
                                    .entries(scope.data);

                                var keyValue = nest.map(function(d) {
                                    return d.key;
                                });

                                if (dimSetting[dim].dimType === 'semiOrdinal') {

                                    keyValue.sort();
                                }

                                dimSetting[dim].keyValue = initializeKeyValueObject(keyValue);


                            };

                            var initializeKeyValueObject = function(keyValue) {

                                var keyObject = {};

                                keyValue.forEach(function(d, i) {
                                    keyObject[d] = {};
                                    keyObject[d].keyValue = d;
                                    keyObject[d].sortedID = i;
                                    keyObject[d].isMinimized = false;
                                    keyObject[d].isMaximized = false;
                                    keyObject[d].calculatedPosition = i;
                                });

                                return keyObject;

                            };


                            var identifyDimDataType = function(dim) {

                                if (isFirstSampleNumber(dim)) {

                                    return identifyOrdinalDimDataType(dim);
                                } else {

                                    return "nominal";
                                }

                            };

                            var identifyOrdinalDimDataType = function(dim) {

                                if (isSemiOrdinalDim(dim)) {

                                    return "semiOrdinal";
                                } else {

                                    return "ordinal";
                                }

                            };

                            var isSemiOrdinalDim = function(dim) {

                                if (getRawNumberOfKeys(dim) < 60) {
                                    return true;
                                } else {
                                    return false;
                                }


                            };

                            var getRawNumberOfKeys = function(dim) {

                                if (!dim) {

                                    return 1;
                                }

                                var nest = d3.nest()
                                    .key(function(d) {
                                        return d[dim];
                                    })
                                    .entries(scope.data);

                                var keyValue = nest.map(function(d) {
                                    return d.key;
                                });

                                return keyValue.length;

                            };

                            var getKeys = function(dim) {

                                if (!dim) {

                                    return [''];
                                }


                                return d3.map(dimSetting[dim].keyValue).keys();
                            };


                            var isFirstSampleNumber = function(dim) {

                                return !isNaN(scope.data[0][dim]);

                            };

                            scope.renderDataChange = function(data, config) {

                                if (!data) {
                                    return;
                                }

                                renderData = data;

                                reloadDataToSVG();

                                identifyAndUpdateDimDataType();

                                scope.handleConfigChange(data, config);

                            }; //End Data change renderer



                            // define render function
                            scope.handleConfigChange = function(data, config) {

                                if (!data) {
                                    return;
                                }

                                renderConfigChange(data, config);

                                handleLensChange(config);

                            };

                            var redrawLensTopic = function(lensInfo) {

                                var itemsOnLens = getLensData(lensInfo);

                                calculatePositionForLensTopic(itemsOnLens, lensInfo);

                                // drawLensItems(itemsOnLens, lensInfo);

                            };

                            var redrawLensRect = function(lensInfo) {

                                var itemsOnLens = getLensData(lensInfo);

                                calculatePositionForLensRect(itemsOnLens, lensInfo);

                                drawLensItems(itemsOnLens, lensInfo);

                            };

                            var redrawLensHist = function(lensInfo) {

                                var itemsOnLens = getLensData(lensInfo);

                                calculatePositionForLensHist(itemsOnLens, lensInfo);

                                drawLensItems(itemsOnLens, lensInfo);

                            };

                            var redrawLensPie = function(lensInfo) {

                                var itemsOnLens = getLensData(lensInfo);

                                calculatePositionForLensPie(itemsOnLens, lensInfo);

                                drawLensItems(itemsOnLens, lensInfo);

                            };



                            var rescale = function(mappedX, width, height) {
                                var maxX = -1e10,
                                    maxY = -1e10,
                                    minX = 1e10,
                                    minY = 1e10,
                                    displayRatio = 0.8;

                                for (var i = 0; i < mappedX.length; i++) {
                                    if (maxX < mappedX[i][0]) maxX = mappedX[i][0];
                                    else if (minX > mappedX[i][0]) minX = mappedX[i][0];
                                    if (maxY < mappedX[i][1]) maxY = mappedX[i][1];
                                    else if (minY > mappedX[i][1]) minY = mappedX[i][1];
                                }

                                var mappedX_tmp = [],
                                    rangeX = maxX - minX,
                                    rangeY = maxY - minY;
                                var rescaleX = (width/rangeX) * displayRatio,
                                    rescaleY = (height/rangeY) * displayRatio;

                                for (var i=0; i<mappedX.length; i++) {
                                    mappedX_tmp.push([mappedX[i][0]*rescaleX, mappedX[i][1]*rescaleY]);
                                }

                                return mappedX_tmp;
                            }


                            var mixTopicColor = function(items,main_k,sub_k,subCluster) {
                                var originalTopicNum = [];
                                var originalTopicColor = [];
                                var subTopicColor = [];

                                for(var i=0;i<main_k;i++) {
                                    var mainColor = color(i);
                                    originalTopicColor.push([   // r,g,b to integer
                                        parseInt(mainColor.slice(1,3), 16),
                                        parseInt(mainColor.slice(3,5), 16), 
                                        parseInt(mainColor.slice(5,7), 16)
                                    ]);
                                }
                                for(var i=0;i<sub_k;i++) {
                                    var originalTopicNum_sub = [];
                                    for(var j=0;j<main_k;j++)  {
                                        originalTopicNum_sub.push(0);
                                    }
                                    originalTopicNum.push(originalTopicNum_sub);
                                }
                                for(var i=0;i<items.length;i++) {
                                    originalTopicNum[items[i].subtopic][parseInt(items[i].cluster)-1]+=1;
                                }
                                for(var i=0;i<sub_k;i++) {
                                    var r=0,g=0,b=0;

                                    for(var j=0;j<main_k;j++) {
                                        r += originalTopicNum[i][j]*originalTopicColor[j][0];
                                        g += originalTopicNum[i][j]*originalTopicColor[j][1];
                                        b += originalTopicNum[i][j]*originalTopicColor[j][2];
                                    }

                                    var original_ith_topic_num = originalTopicNum[i].reduce(function(a,b) { return a+b; });
                                    var range = 100;
                                    r = Math.floor(r/original_ith_topic_num+(Math.random()*range-range/2));
                                    g = Math.floor(g/original_ith_topic_num+(Math.random()*range-range/2));
                                    b = Math.floor(b/original_ith_topic_num+(Math.random()*range-range/2));

                                    r = r>255 ? 255:r;
                                    r = r<0 ? 0:r;
                                    g = g>255 ? 255:g;
                                    g = g<0 ? 0:g;
                                    b = b>255 ? 255:b;
                                    b = b<0 ? 0:b;

                                    r = r<16 ? '0'+r.toString(16):r.toString(16);
                                    g = g<16 ? '0'+g.toString(16):g.toString(16);
                                    b = b<16 ? '0'+b.toString(16):b.toString(16);
                                    subTopicColor.push('#'+r+g+b);
                                }
                                for(var i=0;i<items.length;i++) items[i].subColor = subTopicColor[items[i].subtopic];
                                for(var i=0;i<sub_k;i++) subCluster[i].subColor = subTopicColor[i];
                            }


                            var calculatePositionForLensTopic = function(items, lensInfo) {
                                Y = undefined;

                                var socketId = Math.random().toString(16)
                                var selectedItems = items.map(function(d) {
                                    return +d.filenumber;
                                });

                                socket.on('result data'+socketId, function(data) {
                                        clearInterval(tsne_animation);
                                        var cl_idx_sub = [];
                                        for(var i=0;i<data.cl_idx_sub.length;i++) cl_idx_sub.push(data.cl_idx_sub[i]-1);

                                        for(var i=0;i<items.length;i++) items[i].subtopic=cl_idx_sub[i];
                                        var sub_k = data.Wtopk_sub.length;

                                        subCluster = new Array(sub_k);
                                        for(var i=0;i<subCluster.length;i++) { 
                                            subCluster[i] = {};
                                            subCluster[i].keywords = data.Wtopk_sub[i].slice(0,3);
                                        }
                                        
                                        mixTopicColor(items,10,sub_k,subCluster);
                                        
                                        var distanceMatrix_sub = data.distanceMatrix
                                        var coord = [];
                                        var totalData = scope.data;
                                        for(var i=0;i<totalData.length;i++) {
                                            coord.push([parseFloat(totalData[i].X), parseFloat(totalData[i].Y)]);
                                        }

                                        var N = selectedItems.length;
                                        if (N != distanceMatrix_sub.length){
                                            console.log("error");
                                            return ;    
                                        }
                                        // average distance of the matrix                                        
                                        var totalsum = 0;
                                        for (var i=0; i<N; i++){
                                            for (var j=0; j<N; j++){
                                                totalsum += distanceMatrix_sub[i][j];
                                            }
                                        }
                                        var avg=totalsum/(N*N/2);
                                        
                                        //calculate the centroid coordinate
                                        var ctrary = [];
                                        for(var i=0;i<sub_k;i++) {
                                            ctrary.push([0,0]);
                                        }
                                        var topicNum = new Array(sub_k);
                                        for(var i=0;i<sub_k;i++) topicNum[i]=0;

                                        var summ1 = 0,
                                            summ2 = 0;
                                        for(var i=0;i<coord.length;i++) {
                                            summ1 += coord[i][0];
                                            summ2 += coord[i][1];
                                        }
                                        for(var i=0;i<coord.length;i++) {
                                            coord[i][0] = coord[i][0] - summ1/coord.length;
                                            coord[i][1] = coord[i][1] - summ2/coord.length;
                                        }

                                        for(var i=0;i<selectedItems.length;i++) {
                                            var topicIndex = cl_idx_sub[i];
                                            ctrary[topicIndex][0] += coord[selectedItems[i]-1][0];
                                            ctrary[topicIndex][1] += coord[selectedItems[i]-1][1];
                                            topicNum[topicIndex] += 1;
                                        }
                                        for(var i=0;i<sub_k;i++) {
                                            ctrary[i][0]/=topicNum[i];
                                            ctrary[i][1]/=topicNum[i];
                                        }
                                        var sum1 = 0,
                                            sum2 = 0;
                                        for(var i=0; i<sub_k;i++){
                                            sum1 += ctrary[i][0]
                                            sum2 += ctrary[i][1]
                                        }
                                        var mean1 = sum1/sub_k,
                                            mean2 = sum2/sub_k;

                                        var var1 = 0,
                                            var2 = 0;
                                        for(var i=0; i<sub_k;i++){
                                            var1 += Math.pow((ctrary[i][0]-mean1),2);
                                            var2 += Math.pow((ctrary[i][1]-mean2),2);
                                        }
                                        var1 = var1/sub_k;
                                        var2 = var2/sub_k;
                                        var std1 = Math.sqrt(var1),
                                            std2 = Math.sqrt(var2);

                                        for(var i=0; i<sub_k;i++){
                                            ctrary[i][0] = (ctrary[i][0]-mean1)/std1;
                                            ctrary[i][1] = -(ctrary[i][1]-mean2)/std2;
                                        }

                                        var LM = Math.floor(0.6*N);
                                        if (Y == undefined){
                                            tsne.initDataDist(distanceMatrix_sub,avg, LM);                                            
                                        } else {
                                            tsne.noninitDataDist(distanceMatrix_sub,avg,Y, LM);
                                        }

                                        for (var i = 0; i < 100; i++) tsne.step(sub_k,cl_idx_sub,ctrary, LM);
                                        var intervalNum = 200;
                                        tsne_animation = setInterval(function() {
                                            for (var i = 0; i < 10; i++) tsne.step(sub_k,cl_idx_sub,ctrary, LM);
                                            Y = rescale(tsne.getSolution(), lensInfo.width, lensInfo.height);

                                            calculatePositionUsingSubClusterForLensTopic(items, lensInfo);
                                            drawLensItems(items, lensInfo);

                                            for(var i=0;i<subCluster.length;i++) { 
                                                subCluster[i].X = 0;
                                                subCluster[i].Y = 0;
                                                subCluster[i].num = 0;
                                                subCluster[i].id = i;
                                            }
                                            
                                            for(var i=0;i<items.length;i++) {
                                                var clusterIndex = parseInt(items[i].subtopic);
                                                subCluster[clusterIndex].X += items[i].lensX;
                                                subCluster[clusterIndex].Y += items[i].lensY;
                                                subCluster[clusterIndex].num += 1;
                                            }
                                            
                                            var leftNum = 0;
                                            var gapBetweenLens = 15;
                                            for(var i=0;i<subCluster.length;i++) {
                                                var num = subCluster[i].num;
                                                subCluster[i].X/=num;
                                                subCluster[i].Y/=num;

                                                if(subCluster[i].X<lensInfo.centerX) {
                                                    subCluster[i].textX = lensInfo.centerX-lensInfo.width/2 - gapBetweenLens;
                                                    leftNum+=1
                                                } else {
                                                    subCluster[i].textX = lensInfo.centerX+lensInfo.width/2 + gapBetweenLens;
                                                }
                                            }

                                            var leftGap = lensInfo.height/leftNum,
                                                rightGap = lensInfo.height/(subCluster.length-leftNum);
                                            var leftIndex=1, rightIndex=1;

                                            subCluster.sort(function(a,b) { return a.Y<b.Y ? 1:-1; });

                                            for(var i=0;i<subCluster.length;i++) {
                                                if(subCluster[i].X<lensInfo.centerX) {
                                                    subCluster[i].textY = lensInfo.centerY+lensInfo.height/2+leftGap/2-leftGap*leftIndex;
                                                    leftIndex+=1; 
                                                } else {
                                                    subCluster[i].textY = lensInfo.centerY+lensInfo.height/2+rightGap/2-rightGap*rightIndex;
                                                    rightIndex+=1;
                                                }
                                            }

                                            subCluster.sort(function(a,b) { return a.id>b.id ? 1:-1; });

                                            nodeGroup.select("g.subTopic").remove();
                                            var subTopic = nodeGroup.append("g").attr("class","subTopic");

                                            subTopic.selectAll("text").remove();
                                            subTopic.selectAll("text")
                                                .data(subCluster).enter()
                                                .append("text")
                                                .attr('x', function(d) { return d.textX; })
                                                .attr('y', function(d) { return d.textY; })
                                                .style("fill", function(d) { return d.subColor; })
                                                .text(function(d) {
                                                    return d.keywords.join(' ');
                                                })
                                                .attr('font-family','sans-serif')
                                                .attr('text-anchor', function(d) {
                                                    if(d.X<lensInfo.centerX) return 'end';
                                                    else return 'start'
                                                })
                                                .style('font-size', '10px')
                                                .style('font-weight', 'bold')
                                                .style('pointer-events', 'none');


                                            subTopic.selectAll("line").remove();
                                            subTopic.selectAll("line")
                                                .data(subCluster).enter()
                                                .append("line")
                                                .style("stroke", "black")
                                                .style('stroke-width',0.5)
                                                .attr("x1", function(d) { return d.textX; })
                                                .attr("y1", function(d) { return d.textY; })
                                                .attr("x2", function(d) { return d.X; })
                                                .attr("y2", function(d) { return d.Y; });

                                            var subTexts = subTopic.selectAll("text")[0];
                                            var subTextNum = subTexts.length;
                                            var bboxList = [];
                                            for(var i=0;i<subTextNum;i++) {
                                                bboxList.push(subTexts[i].getBBox());
                                            }

                                            var margin = 7;
                                            subTopic.selectAll("rect").remove();
                                            subTopic.selectAll("rect")
                                                .data(bboxList).enter()
                                                .append("rect")
                                                .attr("x", function(d) { return d.x-margin/2; })
                                                .attr("y", function(d) { return d.y-margin/2; })
                                                .attr("width", function(d) { return d.width+margin; })
                                                .attr("height", function(d) { return d.height+margin; })
                                                .attr("fill", "white")
                                                .style("stroke", "black");


                                            d3.selection.prototype.moveToFront = function() {  
                                                return this.each(function() { 
                                                    this.parentNode.appendChild(this); 
                                                });
                                            };

                                            subTexts.forEach(function(d) {  
                                                d3.select(d).moveToFront();
                                            })

                                            intervalNum-=1;
                                            if(intervalNum==0) clearInterval(tsne_animation);
                                        }, 50);
                                    });

                                    socket.emit('get_subTopic',{'idx':selectedItems, 'socketId':socketId});
                            };


                            var calculatePositionUsingSubClusterForLensTopic = function(items, lensInfo) {
                                var nestedLensItems = d3.nest()
                                    .key(function(d) {
                                        return d.subtopic;
                                    })
                                    .sortKeys(d3.ascending)
                                    .sortValues(function(a, b) {
                                        return a[scope.xdim] - b[scope.xdim];
                                    })
                                    .entries(items);

                                
                                nestedLensItems.forEach(function(d) {
                                    handleOffsetHistLens2(d.values);
                                });

                                
                                for (var i = 0; i < items.length; i++) {
                                    items[i].lensX = lensInfo.centerX + Y[i][0];
                                    items[i].lensY = lensInfo.centerY + Y[i][1];
                                }
                            };



                            var calculatePositionForLensRect = function(items, lensInfo) {

                                items.sort(sortFuncByColorDimension());

                                items.forEach(function(d, i) {

                                    d.clusterID = i;
                                    d.lensX = lensInfo.centerX;
                                    d.lensY = lensInfo.centerY;

                                })

                                var box = getLensClusterSize(1, lensInfo);

                                var maxNumElementInCluster = items.length;

                                var size = calculateNodesSizeForAbsolute(box, maxNumElementInCluster);

                                if (size > maxDotSize) {

                                    size = maxDotSize;
                                }

                                handleOffsetTopicLens(items, box, size);

                            };

                            var handleOffsetTopicLens = function(cluster, box, size) {

                                if (box.widthOfBox > box.heightOfBox) {

                                    handleOffsetRectLensHorizontally(cluster, box, size);

                                } else {

                                    handleOffsetRectLensVertically(cluster, box, size);
                                }

                            };

                            var handleOffsetRectLens = function(cluster, box, size) {

                                if (box.widthOfBox > box.heightOfBox) {

                                    handleOffsetRectLensHorizontally(cluster, box, size);

                                } else {

                                    handleOffsetRectLensVertically(cluster, box, size);
                                }

                            };

                            var handleOffsetHistLens = function(cluster, box, size) {

                                if (box.widthOfBox > box.heightOfBox) {

                                    handleOffsetHistLensHorizontally(cluster, box, size);

                                } else {

                                    handleOffsetHistLensVertically(cluster, box, size);
                                }

                            };


                            var handleOffsetHistLens2 = function(cluster) {
                                var size = 3;
                                cluster.forEach(function(d, i, j) {

                                    d.nodeWidthLens = size;
                                    d.nodeHeightLens = size;

                                    d.YOffsetLens = 0;
                                    d.XOffsetLens = 0;

                                });

                            };

                            var calculatePositionForLensHist = function(items, lensInfo) {

                                var nestedLensItems = d3.nest()
                                    .key(function(d) {
                                        return d[scope.config.colorDim];
                                    })
                                    .sortKeys(d3.ascending)
                                    .sortValues(function(a, b) {
                                        return a[scope.xdim] - b[scope.xdim];
                                    })
                                    .entries(items);


                                assignClusterIDOfNodesInOneKeyNestedItems(nestedLensItems);

                                var box = getLensClusterSize(nestedLensItems.length, lensInfo);

                                var maxNumElementInCluster = getClusterWithMaximumPopulationFromOneKeyNestedItems(nestedLensItems);

                                var size = calculateNodesSizeForAbsolute(box, maxNumElementInCluster);

                                if (size > maxDotSize) {

                                    size = maxDotSize;
                                }

                                nestedLensItems.forEach(function(d) {
                                    handleOffsetHistLens(d.values, box, size);
                                });

                                nestedLensItems.forEach(function(d, i) {

                                    d.values.forEach(function(d) {
                                        d.lensX = lensInfo.centerX - lensInfo.width / 2 + (0.5 + i) * box.widthOfBox;
                                        d.lensY = lensInfo.centerY;
                                    });
                                });

                            };

                            var calculatePositionForLensPie = function(items, lensInfo) {

                                items.forEach(function(d, i) {
                                    d.lensX = lensInfo.centerX;
                                    d.lensY = lensInfo.centerY;

                                });

                                var nestedLensItems = d3.nest()
                                    .key(function(d) {
                                        return d[scope.config.colorDim];
                                    })
                                    .sortKeys(d3.ascending)
                                    .sortValues(function(a, b) {
                                        return a[scope.xdim] - b[scope.xdim];
                                    })
                                    .entries(items);

                                var numElement = items.length;

                                var layerSetting = calculateLayerSettingForPieLens(lensInfo, numElement);

                                if (layerSetting.dotR > maxDotSize / 2) {

                                    layerSetting.dotR = maxDotSize / 2;
                                }

                                var clusterAngle = getClusterAngle(nestedLensItems, layerSetting, numElement);

                                nestedLensItems.forEach(function(d, i) {
                                    handleOffsetPieLens(d.values, layerSetting, clusterAngle[i], lensInfo);
                                });

                            };

                            var handleOffsetPieLens = function(items, layerSetting, clusterAngle, lensInfo) {

                                var currentLayer = 0;

                                var angleOffset = clusterAngle.startAngle;

                                items.forEach(function(d, i, j) {

                                    d.nodeWidthLens = layerSetting.dotR * 2;
                                    d.nodeHeightLens = layerSetting.dotR * 2;

                                    angleOffset = angleOffset + layerSetting.layerIncrementalAngle[currentLayer];

                                    if (angleOffset >= clusterAngle.endAngle) {

                                        angleOffset = clusterAngle.startAngle;
                                        currentLayer++;
                                    }

                                    var angle = angleOffset;
                                    var innerR = layerSetting.layerInnerRadius[currentLayer];

                                    d.XOffsetLens = innerR * Math.cos(angle);
                                    d.YOffsetLens = -1 * innerR * Math.sin(angle);


                                });
                            }

                            var getClusterAngle = function(nestedItems, layerSetting, numElement) {

                                var clusterNumber = nestedItems.length;
                                var offsetAngle = 0;

                                var clusterAngle = nestedItems.map(function(d, i, j) {

                                    var startAngle = offsetAngle;
                                    var endAngle = startAngle + 2 * Math.PI * (d.values.length / numElement);
                                    offsetAngle = endAngle;

                                    return {
                                        'startAngle': startAngle,
                                        'endAngle': endAngle
                                    };
                                });

                                return clusterAngle;

                            };

                            var calculateLayerSettingForPieLens = function(lensInfo, numElement) {

                                var numLayer = 1;

                                while (!isNumLayerLargeEnough(numLayer, lensInfo, numElement)) {
                                    numLayer++;
                                }

                                return calculateLayerSettingForPieLensWithNumLayer(numLayer, lensInfo, numElement);
                            };

                            var calculateLayerSettingForPieLensWithNumLayer = function(numLayer, lensInfo, numElement) {

                                var layerSetting = {};

                                layerSetting.numLayer = numLayer;
                                layerSetting.dotR = getDotRadiusFromNumLayer(numLayer, lensInfo);

                                layerSetting.layerInnerRadius = [];
                                layerSetting.layerIncrementalAngle = [];
                                layerSetting.itemSetting = [];
                                var itemCountStart = 0;

                                for (var layer = 1; layer <= layerSetting.numLayer; layer++) {

                                    var innerR = getInnerRadiusPieLens(lensInfo, layer, layerSetting.dotR);
                                    layerSetting.layerInnerRadius.push(innerR);

                                    var incrementalAngle = getIncrementalAngleForPieLens(layerSetting.dotR, lensInfo, layer);
                                    layerSetting.layerIncrementalAngle.push(incrementalAngle);

                                    // for (var itemCount = itemCountStart; itemCount <= itemCountStart + count; itemCount++) {

                                    //     var itemSetting = {}

                                    //     itemSetting.layer = layer;
                                    //     itemSetting.incrementalAngle = incrementalAngle;
                                    //     layerSetting.itemSetting[itemCount] = itemSetting;

                                    //     console.log(itemCount);

                                    // }

                                    // itemCountStart = count; 
                                }

                                for (var itemCount = 0; itemCount < numElement; itemCount++) {


                                }

                                return layerSetting;

                            };

                            var isNumLayerLargeEnough = function(numLayer, lensInfo, numElement) {

                                var dotR = getDotRadiusFromNumLayer(numLayer, lensInfo);

                                if (dotR >= lensInfo.innerRadius) {

                                    return false;
                                }

                                var accumulatedItemsCount = 0;

                                for (var i = 1; i <= numLayer; i++) {

                                    accumulatedItemsCount = accumulatedItemsCount + numItemsInSingleLayer(i, dotR, lensInfo);
                                }

                                return (numElement <= accumulatedItemsCount);
                            };

                            var numItemsInSingleLayer = function(layerCount, dotR, lensInfo) {

                                var angleForDot = getIncrementalAngleForPieLens(dotR, lensInfo, layerCount);
                                var numItems = Math.floor(2 * Math.PI / angleForDot);

                                return numItems;
                            };

                            var getIncrementalAngleForPieLens = function(dotR, lensInfo, layerCount) {

                                var innerRadius = getInnerRadiusPieLens(lensInfo, layerCount, dotR);

                                var halfAngle = Math.PI / 2 - Math.acos(dotR / innerRadius);

                                return halfAngle * 2;

                            };

                            var getInnerRadiusPieLens = function(lensInfo, layerCount, dotR) {

                                return lensInfo.innerRadius + (2 * (layerCount - 1) * dotR);
                            }

                            var getDotRadiusFromNumLayer = function(numLayer, lensInfo) {

                                var dotR = (lensInfo.outerRadius - lensInfo.innerRadius) / (2 * numLayer - 1);

                                return dotR;
                            };

                            var getLensClusterSize = function(clusterNumber, lensInfo) {

                                var lensClusterSize = {};

                                lensClusterSize.widthOfBox = lensInfo.width / clusterNumber;

                                lensClusterSize.heightOfBox = lensInfo.height;

                                return lensClusterSize;
                            }

                            var drawLensItems = function(itemsOnLens, lensInfo) {
                                nodeGroup.selectAll(".dot")
                                    .data(itemsOnLens, function(d) {
                                        return d.id;
                                    }).remove();

                                var lensItems = nodeGroup.selectAll(".lensItems")
                                    .data(itemsOnLens, function(d) {
                                        return d.id;
                                    });


                                //Update
                                //Transition from previous place to new place
                                lensItems.transition()
                                    .attr("width", function(d) {
                                        return +d.nodeWidthLens*1.5;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeightLens*1.5;
                                    })
                                    .attr("x", function(d) {
                                        return d.lensX;
                                    })
                                    .attr("y", function(d) {
                                        return d.lensY;
                                    })
                                    .attr('fill-opacity', 0.8)
                                    .style("fill", function(d) {
                                        return d.subColor;
                                    })
                                .attr("transform", function(d, i) {
                                    return "translate(" + (d.XOffsetLens) + "," + (-(d.YOffsetLens)) + ") ";
                                });


                                //Enter
                                //Append new circle
                                //Transition from Original place to new place
                                lensItems.enter().append("rect")
                                    .classed({
                                        'lensItems': true,
                                        'dot': false
                                    })
                                    .attr("y", function(d,i) { 
                                        return yMap(d); })
                                    .attr("x", function(d) {
                                        return xMap(d); })
                                    .attr("width", function(d) {
                                        return +d.nodeWidth*1.5;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeight*1.5;
                                    })
                                    .attr("rx", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .attr("ry", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .style("fill", function(d) {
                                        return color(d.subtopic+10);
                                    })
                                    .transition()
                                    .attr("x", function(d) {
                                        return d.lensX;
                                    })
                                    .attr("y", function(d) {
                                        return d.lensY;
                                    })
                                    .attr("width", function(d) {
                                        return +d.nodeWidthLens;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeightLens;
                                    })
                                    .attr("transform", function(d, i) {
                                        return "translate(" + (d.XOffsetLens) + "," + (-(d.YOffsetLens)) + ") ";
                                    });


                                lensItems.on("mouseover", function(d) {
                                    var pageX = d3.event.pageX,
                                        pageY = d3.event.pageY;

                                    mouseOverTimer = setTimeout(function() {
                                        d3.select('g.subTopic').selectAll('rect').style('opacity',0);
                                        d3.select('g.subTopic').selectAll('line').style('opacity',0);
                                        d3.select('g.subTopic').selectAll('text').style('opacity',0);

                                        var authorList = d.authors.split(', ');
                                        var tooltipString = '<b>Authors:<br/>';
                                        for(var i=0;i<authorList.length;i++) {
                                            if(i==4) tooltipString+='<br/>'+authorList[i];
                                            else tooltipString+=authorList[i];
                                        }
                                        tooltipString+='<br/>Title:<br/>';
                                        var titleWords = d.title.split(' ');
                                        var accumulatedLength=0;
                                        var maxLineLength = 40;
                                        for(var i=0;i<titleWords.length;i++) {
                                            accumulatedLength+=titleWords[0].length+1;
                                            if(accumulatedLength>maxLineLength) {
                                                accumulatedLength=0;
                                                tooltipString+='<br/>'+titleWords[i]+' ';
                                            }
                                            else tooltipString+=titleWords[i]+' ';
                                        }
                                        tooltipString+='</b>'

                                        tooltip.transition()
                                            .duration(100)
                                            .style("opacity", 0.9);

                                        tooltip.html(tooltipString)
                                            .style("left", (pageX + 5) + "px")
                                            .style("top", (pageY - 28) + "px");
                                    }, 700);

                                })
                                .on("mouseout", function(d) {
                                    clearTimeout(mouseOverTimer);
                                    d3.select('g.subTopic').selectAll('rect').style('opacity',1);
                                    d3.select('g.subTopic').selectAll('line').style('opacity',1);
                                    d3.select('g.subTopic').selectAll('text').style('opacity',1);
                                    tooltip.style("opacity", 0);
                                });

                            
                                //Exit
                                //Transition from previous place to original place
                                //remove circle
                                lensItems.exit()
                                    .classed({
                                        'dot': true,
                                        'lensItems': false
                                    })
                                    .transition()
                                    .duration(300)
                                    .attr("width", function(d) {
                                        return +d.nodeWidth*1.5;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeight*1.5;
                                    })
                                    .attr("y", function(d) {
                                        return yMap(d); })
                                    .attr("x", function(d) {
                                        return xMap(d); })
                                    .attr("transform", function(d, i) {
                                        return "translate(" + (d.XOffset) + "," + (-(d.YOffset)) + ") ";
                                    })
                                    .style('fill', function(d) {
                                        return color(d[scope.config.colorDim]-1);
                                    });
                            };

                            var reverseTransform = function(coor, transfactor, scalefactor) {
                                return coor;
                            };



                            var getLensData = function(lensInfo) {

                                var itemsOnLens = scope.data.filter(function(d) {

                                    if (xMap(d) < reverseTransform(lensInfo.centerX + lensInfo.width / 2, scope.context.translate[0], scope.context.scale) && xMap(d) > reverseTransform(lensInfo.centerX - lensInfo.width / 2, scope.context.translate[0], scope.context.scale)) {

                                        if (yMap(d) < reverseTransform(lensInfo.centerY + lensInfo.height / 2, scope.context.translate[1], scope.context.scale) && yMap(d) > reverseTransform(lensInfo.centerY - lensInfo.height / 2, scope.context.translate[1], scope.context.scale)) {

                                            return d;
                                        }
                                    }

                                });

                                if (lensInfo.type === 'pie') {
                                    itemsOnLens = itemsOnLens.filter(function(d) {

                                        var x = xMap(d) - lensInfo.centerX;
                                        var y = yMap(d) - lensInfo.centerY;

                                        var dist = Math.sqrt(x * x + y * y);
                                        if (dist < lensInfo.outerRadius) {
                                            return d;
                                        }
                                    });
                                }

                                //console.log(itemsOnLens);


                                return itemsOnLens;


                            };

                            var drawInitialTopicLensItems = function(centerX, centerY, width, height) {

                                //var lensInfo = {};

                                lensInfo.centerX = centerX;
                                lensInfo.centerY = centerY;
                                lensInfo.type = 'topic';
                                lensInfo.width = width;
                                lensInfo.height = height;
                                redrawLensTopic(lensInfo);
                            };

                            var drawInitialRectLensItems = function(centerX, centerY, width, height) {

                                var lensInfo = {};

                                lensInfo.centerX = centerX;
                                lensInfo.centerY = centerY;
                                lensInfo.type = 'rect';
                                lensInfo.width = width;
                                lensInfo.height = height;
                                redrawLensRect(lensInfo);
                            };

                            var drawInitialHistLensItems = function(centerX, centerY, width, height) {

                                var lensInfo = {};

                                lensInfo.centerX = centerX;
                                lensInfo.centerY = centerY;
                                lensInfo.type = 'hist';
                                lensInfo.width = width;
                                lensInfo.height = height;

                                redrawLensHist(lensInfo);
                            };

                            var drawInitialPieLensItems = function(centerX, centerY, width, height) {

                                var lensInfo = {};

                                lensInfo.centerX = centerX;
                                lensInfo.centerY = centerY;
                                lensInfo.type = 'pie';
                                lensInfo.width = width;
                                lensInfo.height = height;
                                lensInfo.outerRadius = initialLensSize / 2;
                                lensInfo.innerRadius = initialInnerRadiusOfPieLens;


                                redrawLensPie(lensInfo);
                            };

                            var dragmoveTopicLens = function() {

                                var xPos, yPos;

                                // var lensInfo = {};

                                d3.select(this)
                                    .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.x)) - initialHistLensWidth / 2)
                                    .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.y)) - initialHistLensHeight / 2);

                                // labelDiv.text(xPos);

                                lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                lensInfo.type = 'hist';
                                lensInfo.width = initialHistLensWidth;
                                lensInfo.height = initialHistLensHeight;


                                // redrawLensTopic(lensInfo);

                            };

                            var dragendTopicLens = function() {

                                var xPos, yPos;

                                // var lensInfo = {};

                                //d3.select(this)
                                //    .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.sourceEvent.x)) - initialHistLensWidth / 2)
                                //    .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.sourceEvent.y)) - initialHistLensHeight / 2);

                                // // labelDiv.text(xPos);

                                //lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                //lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                //lensInfo.type = 'hist';
                                //lensInfo.width = initialHistLensWidth;
                                //lensInfo.height = initialHistLensHeight;


                                redrawLensTopic(lensInfo);

                            };


                            var dragmoveRectLens = function() {

                                var xPos, yPos;

                                var lensInfo = {};

                                d3.select(this)
                                    .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.x)) - initialHistLensWidth / 2)
                                    .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.y)) - initialHistLensHeight / 2);

                                // labelDiv.text(xPos);

                                lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                lensInfo.type = 'hist';
                                lensInfo.width = initialHistLensWidth;
                                lensInfo.height = initialHistLensHeight;


                                // redrawLensTopic(lensInfo);

                            };

                            var dragmoveHistLens = function() {

                                var xPos, yPos;

                                // var lensInfo = {};

                                d3.select(this)
                                    .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.x)) - initialHistLensWidth / 2)
                                    .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.y)) - initialHistLensHeight / 2);

                                // labelDiv.text(xPos);

                                lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                lensInfo.type = 'hist';
                                lensInfo.width = initialHistLensWidth;
                                lensInfo.height = initialHistLensHeight;


                                redrawLensHist(lensInfo);

                            };

                            var dragendHistLens = function() {

                                // var xPos, yPos;

                                // var lensInfo = {};

                                // d3.select(this)
                                //     .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.x)) - initialHistLensWidth / 2)
                                //     .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.y)) - initialHistLensHeight / 2);

                                // // labelDiv.text(xPos);

                                // lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                // lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                // lensInfo.type = 'hist';
                                // lensInfo.width = initialHistLensWidth;
                                // lensInfo.height = initialHistLensHeight;


                                redrawLensHist(lensInfo);

                            };

                            var dragmovePieLens = function() {

                                var xPos, yPos;

                                var lensInfo = {};

                                d3.select(this)
                                    .attr("x", xPos = Math.max(initialLensSize / 2, Math.min(width - initialLensSize / 2, d3.event.x)) - initialLensSize / 2)
                                    .attr("y", yPos = Math.max(initialLensSize / 2, Math.min(height - initialLensSize / 2, d3.event.y)) - initialLensSize / 2);

                                // labelDiv.text(xPos);

                                lensInfo.centerX = xPos + initialLensSize / 2;
                                lensInfo.centerY = yPos + initialLensSize / 2;
                                lensInfo.type = 'pie';
                                lensInfo.width = initialLensSize;
                                lensInfo.height = initialLensSize;

                                lensInfo.outerRadius = initialLensSize / 2;
                                lensInfo.innerRadius = initialInnerRadiusOfPieLens;

                                redrawLensPie(lensInfo);


                            };

                            var dragmoveCircle = function() {

                                var xPos, yPos;

                                d3.select(this)
                                    .attr("cx", xPos = Math.max(initialLensSize, Math.min(width - initialLensSize, d3.event.x)))
                                    .attr("cy", yPos = Math.max(initialLensSize, Math.min(height - initialLensSize, d3.event.y)));

                                // labelDiv.text(xPos);

                            }

                            var lensInfo = {};

                            var handleTopicLensChange = function() {

                                clearLens();

                                //var lensInfo = {};

                                  var drag = d3.behavior.drag()
                                    .on("drag", dragmoveTopicLens)
                                    .on("dragstart", function() {
                                        d3.event.sourceEvent.stopPropagation(); // silence other listeners
                                        clearInterval(tsne_animation);
                                        tsne_animation = null;
                                    })
                                    .on("dragend", dragendTopicLens);


                                nodeGroup.append("rect")
                                    .attr("class", "lens")
                                    .attr("x", width / 2)
                                    .attr("y", height / 2)
                                    .attr("width", initialHistLensWidth)
                                    .attr("height", initialHistLensHeight)
                                    .call(drag)
                                
                                d3.select('body').on('keydown',function(){
                                    if(d3.event.keyCode==16) {
                                        shiftDowned=true;
                                    }
                                    else if(shiftDowned&&d3.event.keyCode==38){
                                        var lens = nodeGroup.select('.lens');
                                        var width = parseFloat(lens.attr('width')),
                                            height = parseFloat(lens.attr('height')),
                                            newWidth = width*1.2,
                                            newHeight = height*1.2;
                                        lens.attr('width',function(){ return newWidth; })
                                            .attr('height',function(){ return newHeight; });
                                        initialHistLensWidth = newWidth;
                                        initialHistLensHeight = newHeight;
                                        lensInfo.centerX += width*0.1;
                                        lensInfo.centerY += height*0.1;
                                        lensInfo.width = newWidth;
                                        lensInfo.height = newHeight;
                                        clearInterval(tsne_animation);
                                        redrawLensTopic(lensInfo);
                                    }
                                    else if(shiftDowned&&d3.event.keyCode==40){
                                        var lens = nodeGroup.select('.lens')
                                        var width = parseFloat(lens.attr('width')),
                                            height = parseFloat(lens.attr('height')),
                                            newWidth = width*0.8,
                                            newHeight = height*0.8;
                                        lens.attr('width',function(){ return newWidth; })
                                            .attr('height',function(){ return newHeight; });
                                        initialHistLensWidth = newWidth;
                                        initialHistLensHeight = newHeight;
                                        lensInfo.centerX -= width*0.1;
                                        lensInfo.centerY -= height*0.1;
                                        lensInfo.width = newWidth;
                                        lensInfo.height = newHeight;
                                        clearInterval(tsne_animation);
                                        redrawLensTopic(lensInfo);
                                    }
                                })
                                .on('keyup',function(event){
                                    shiftDowned=false;
                                });

                                drawInitialTopicLensItems(width / 2 + initialHistLensWidth / 2, height / 2 + initialHistLensHeight / 2, initialHistLensWidth, initialHistLensHeight);
                            };

                            var handleRectLensChange = function() {

                                clearLens();

                                var drag = d3.behavior.drag()
                                    .on("drag", dragmoveRectLens)
                                    .on("dragstart", function() {
                                        d3.event.sourceEvent.stopPropagation(); // silence other listeners
                                    });



                                nodeGroup.append("rect")
                                    .attr("class", "lens")
                                    .attr("x", width / 2)
                                    .attr("y", height / 2)
                                    .attr("width", initialLensSize)
                                    .attr("height", initialLensSize)
                                    .call(drag);

                                drawInitialRectLensItems(width / 2 + initialLensSize / 2, height / 2 + initialLensSize / 2, initialLensSize, initialLensSize);
                            };

                            var handleHistLensChange = function() {

                                clearLens();

                                var drag = d3.behavior.drag()
                                    .on("drag", dragmoveHistLens)
                                    .on("dragstart", function() {
                                        d3.event.sourceEvent.stopPropagation(); // silence other listeners
                                    })
                                    .on("dragend", dragendHistLens);

                                nodeGroup.append("rect")
                                    .attr("class", "lens")
                                    .attr("x", width / 2)
                                    .attr("y", height / 2)
                                    .attr("width", initialHistLensWidth)
                                    .attr("height", initialHistLensHeight)
                                    .call(drag);

                                drawInitialHistLensItems(width / 2 + initialHistLensWidth / 2, height / 2 + initialHistLensHeight / 2, initialHistLensWidth, initialHistLensHeight);

                            };

                            var handlePieLensChange = function() {

                                clearLens();

                                var drag = d3.behavior.drag()
                                    .on("drag", dragmovePieLens)
                                    .on("dragstart", function() {
                                        d3.event.sourceEvent.stopPropagation(); // silence other listeners
                                    });

                                nodeGroup.append("rect")
                                    .attr("class", "lens")
                                    .attr("x", width / 2)
                                    .attr("y", height / 2)
                                    .attr("width", initialLensSize)
                                    .attr("height", initialLensSize)
                                    .attr("rx", initialLensSize / 2)
                                    .attr("ry", initialLensSize / 2)
                                    .call(drag);

                                drawInitialPieLensItems(width / 2 + initialHistLensWidth / 2, height / 2 + initialHistLensHeight / 2, initialHistLensWidth, initialHistLensHeight);

                            };

                            var handleLensChange = function(config) {

                                if (config.lens === "noLens" || config.isGather !== 'scatter') {

                                    clearLens();

                                } else if (config.lens === "rectLens") {

                                    handleRectLensChange();

                                } else if (config.lens === "histLens") {

                                    handleHistLensChange();

                                } else if (config.lens === "pieLens") {

                                    handlePieLensChange();
                                } else if (config.lens === "topicLens") {

                                    handleTopicLensChange();
                                }

                            };

                            var clearLens = function() {

                                nodeGroup.selectAll(".lens").remove();
                                nodeGroup.selectAll(".lensItems")
                                    .classed({
                                        'dot': true,
                                        'lensItems': false
                                    })
                                    .transition()
                                    .duration(500)
                                    .attr("width", function(d) {
                                        // console.log(initialSquareLenth);
                                        return +d.nodeWidth;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeight;
                                    })
                                    .attr("y", yMap)
                                    .attr("x", xMap)
                                    .attr("transform", function(d, i) {
                                        return "translate(" + (d.XOffset) + "," + (-(d.YOffset)) + ") ";
                                    });
                            }

                            var updateSizeSVG = function(config) {
                                // XPadding = 60;
                                // YPadding = 30;
                                //Update size of SVG

                                if (scope.config.matrixMode === false) {
                                    outerWidth = d3.select(iElement[0]).node().offsetWidth;
                                } else {
                                    outerWidth = d3.select(".matrixGroup").node().offsetWidth;

                                    outerWidth = outerWidth / (scope.config.dims.length) - 2;

                                }
                                // calculate the height
                                outerHeight = outerWidth / config.SVGAspectRatio;

                                svg.attr('height', outerHeight)
                                    .attr('width', outerWidth);

                                width = outerWidth - 2 * margin;
                                height = outerHeight - 2 * margin;

                            };

                            var renderConfigChange = function(data, config) {


                                updateSizeSVG(config);

                                //Call separate render for the rendering

                                drawPlot();

                                configZoomToolbar();

                                configBrushToolbar();

                                configZoom();

                            };

                            var configZoomToolbar = function() {

                                d3.select("#toolbarPanZoom").on("click", setZoomMode);


                                function setZoomMode() {

                                    configZoom();
                                }


                            };

                            var configBrushToolbar = function() {

                                d3.select("#toolbarSelect").on("click", setSelectMode);

                                function setSelectMode() {

                                    configBrush();
                                }

                            };

                            var configBrush = function() {

                                brush = brushGroup.append("g")
                                    .datum(function() {
                                        return {
                                            selected: false,
                                            previouslySelected: false
                                        };
                                    })
                                    .attr("class", "brush")
                                    .call(d3.svg.brush()
                                        .x(xScale)
                                        .y(yScale)
                                        .on("brushstart", function(d) {
                                            node.each(function(d) {

                                                // if (d.Name.indexOf("ciera") > 0) {
                                                //     console.log(d);
                                                // }

                                                d.previouslySelected = d3.event.sourceEvent.shiftKey && d.selected;
                                            });
                                        })
                                        .on("brush", function() {
                                            var extent = d3.event.target.extent();
                                            node.classed("selected", function(d) {

                                                //     return d.selected = d.previouslySelected ^
                                                //         (xScale(extent[0][0]) <= xMap(d) && xMap(d) < xScale(extent[1][0]) && yScale(extent[0][1]) >= yMap(d) && yMap(d) > yScale(extent[1][1]));
                                                // });

                                                var nodeIndex = scope.dimsum.selectionSpace.indexOf(d.id);

                                                if (d.previouslySelected ^ (xScale(extent[0][0]) <= xMap(d) && xMap(d) < xScale(extent[1][0]) && yScale(extent[0][1]) >= yMap(d) && yMap(d) > yScale(extent[1][1]))) {

                                                    if (nodeIndex == -1) {
                                                        scope.dimsum.selectionSpace.push(d.id);
                                                    }
                                                } else {

                                                    if (nodeIndex != -1) {
                                                        scope.dimsum.selectionSpace.splice(nodeIndex, 1);
                                                    }


                                                }

                                            });
                                            scope.$apply();
                                            scope.handleDimsumChange(scope.dimsum);


                                        })
                                        .on("brushend", function() {
                                            d3.event.target.clear();
                                            d3.select(this).call(d3.event.target);
                                        }));



                                d3.select("#toolbarSelect").classed("active", true);
                                d3.select("#toolbarPanZoom").classed("active", false);




                            };

                            var zoomUsingContext = function() {


                                zoom.translate(scope.context.translate);
                                zoom.scale(scope.context.scale);

                                svgGroup.select(".x.axis").call(xAxis);
                                svgGroup.select(".y.axis").call(yAxis);

                                nodeGroup.attr("transform", "translate(" + scope.context.translate + ")scale(" + scope.context.scale + ")");


                                scope.comment = false;

                            };

                            var configZoom = function() {

                                removeBrushMode();

                                function removeBrushMode() {

                                    d3.selectAll(".brush").remove();


                                }

                                zoom = d3.behavior.zoom()
                                    .x(xScale)
                                    .y(yScale)
                                    .scaleExtent([1, 100])
                                    .on("zoom", zoomed);

                                var panning = d3.behavior.zoom()
                                    .scaleExtent([1, 1])
                                    .x(xScale)
                                    .y(yScale)
                                    .on("zoom", zoomed);

                                svgGroup.call(zoom);
                                var dragmoveRectLens = function() {

                                    var xPos, yPos;

                                    var lensInfo = {};

                                    d3.select(this)
                                        .attr("x", xPos = Math.max(initialHistLensWidth / 2, Math.min(width - initialHistLensWidth / 2, d3.event.x)) - initialHistLensWidth / 2)
                                        .attr("y", yPos = Math.max(initialHistLensHeight / 2, Math.min(height - initialHistLensHeight / 2, d3.event.y)) - initialHistLensHeight / 2);

                                    // labelDiv.text(xPos);

                                    lensInfo.centerX = xPos + initialHistLensWidth / 2;
                                    lensInfo.centerY = yPos + initialHistLensHeight / 2;
                                    lensInfo.type = 'hist';
                                    lensInfo.width = initialHistLensWidth;
                                    lensInfo.height = initialHistLensHeight;


                                    // redrawLensTopic(lensInfo);

                                };


                                function zoomed() {

                                    // scope.$apply();

                                    // scope.comment = false;


                                    // zoom.x(xScale).y(yScale);

                                    svgGroup.select(".x.axis").call(xAxis);
                                    svgGroup.select(".y.axis").call(yAxis);

                                    scope.context.translate = d3.event.translate;
                                    scope.context.scale = d3.event.scale;
                                    scope.context.xDomain = zoom.x().domain();
                                    scope.context.yDomain = zoom.y().domain();

                                    var range = getExtentConsideringXY(scope.xdim, scope.ydim);

                                    xRange = range.xRange;
                                    yRange = range.yRange;

                                    var xScaleForNodes = d3.scale.linear().range([0, width]);
                                    xScaleForNodes.domain(xRange);

                                    xMap = function(d) {
                                        return xScaleForNodes(xValue(d));
                                    };

                                    var yScaleForNodes = d3.scale.linear().range([height, 0]);
                                    yScaleForNodes.domain(yRange);
                                    yMap = function(d) {
                                        return yScaleForNodes(yValue(d));
                                    };



                                    nodeGroup.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
                                    //nodeGroup.attr("transform", "translate(" + d3.event.translate + ")");

                                }

                                function zoomReset() {

                                    svgGroup.select(".x.axis").call(xAxis);
                                    svgGroup.select(".y.axis").call(yAxis);

                                }

                                setClipPathForAxes();

                                d3.select("#toolbarReset").on("click", reset);

                                d3.select("#toolbarSelect").classed("active", false);

                                d3.select("#toolbarPanZoom").classed("active", true);

                                function reset() {

                                    resetSelection();

                                    nodeGroup.transition()
                                        .duration(700)
                                        .attr("transform", "translate(0,0) scale(1)");

                                    scope.context.translate = [0, 0];
                                    scope.context.scale = 1;

                                    // scope.$apply();

                                    d3.transition().duration(700).tween("zoom", function() {

                                        var range = getExtentConsideringXY(scope.xdim, scope.ydim);

                                        var xRange = range.xRange;
                                        var yRange = range.yRange;


                                        if (scope.config.isGather === 'gather') {

                                            var typeOfXYDim = findTypeOfXYDim();

                                            if (typeOfXYDim === 'XNomYOrd') {

                                                var yRange = getExtentFromCalculatedPointsForBinnedGather(scope.ydim);

                                            } else if (typeOfXYDim === 'XOrdYNom') {

                                                var xRange = getExtentFromCalculatedPointsForBinnedGather(scope.xdim);

                                            }

                                        }

                                        var ix = d3.interpolate(xScale.domain(), xRange),
                                            iy = d3.interpolate(yScale.domain(), yRange);

                                        return function(t) {
                                            zoom.x(xScale.domain(ix(t))).y(yScale.domain(iy(t)));

                                            zoomReset();
                                        };
                                    });
                                }

                                if (scope.comment === true) {

                                    zoomUsingContext();

                                } else {

                                    scope.context.translate = [0, 0];
                                    scope.context.scale = 1;
                                    scope.context.xDomain = xScale.domain();
                                    scope.context.yDomain = yScale.domain();

                                }

                            };


                            var resetSelection = function() {

                                if (!scope.dimsum) {
                                    return;

                                }

                                scope.dimsum.selectionSpace = [];
                                scope.handleDimsumChange(scope.dimsum);
                                scope.$apply();
                            }

                            var setClipPathForAxes = function() {

                                var clipXAxis = xAxisNodes.append("clipPath")
                                    .attr("id", "clipXAxis")
                                    .append("rect")
                                    .attr("x", 0)
                                    .attr("y", 0)
                                    .attr("width", width)
                                    .attr("height", margin);

                                xAxisNodes.attr("clip-path", "url(#clipXAxis)");


                                var clipYAxis = yAxisNodes.append("clipPath")
                                    .attr("id", "clipYAxis")
                                    .append("rect")
                                    .attr("x", -300)
                                    .attr("y", -40)
                                    .attr("width", 300)
                                    .attr("height", height + 40);

                                yAxisNodes.attr("clip-path", "url(#clipYAxis)");


                            };

                            var drawPlot = function() {

                                // drawBackground();

                                drawNodes();

                                drawMask();

                                drawAxesAndLegends();

                            };

                            var drawBackground = function() {

                                nodeGroup.append("rect")
                                    .attr("class", "overlay")
                                    .attr("width", width)
                                    .attr("height", height);
                            };

                            var drawMask = function() {

                                var clip = maskGroup.append("clipPath")
                                    .attr("id", "clip")
                                    .append("rect")
                                    .attr("x", 0)
                                    .attr("y", 0)
                                    .attr("width", width)
                                    .attr("height", height);

                                maskGroup.attr("clip-path", "url(#clip)");


                            }

                            var drawAxesAndLegends = function() {

                                if (scope.config.matrixMode === false) {

                                    drawAxes();

                                    // drawLegends();

                                } else {

                                    // drawAxes();

                                    drawBoundaryForMatrix();
                                }
                            }


                            var drawBoundaryForMatrix = function() {

                                svgGroup.selectAll(".matrixFrame").remove();

                                svgGroup.append("rect")
                                    .attr("class", "matrixFrame")
                                    .attr("x", -margin)
                                    .attr("y", -margin)
                                    .attr("width", width + 2 * margin - 2)
                                    .attr("height", height + 2 * margin - 2);


                            };



                            var drawNodesForSameOrdDimGather = function() {

                                prepareScaleForSameOrdDimGather();

                                calculateParametersOfNodesForSameOrdDimGather();

                                drawNodesInSVGForSameOrdDimGather();

                            };

                            var isSameOrdDimGather = function() {

                                if (scope.config.isGather === 'gather' &&
                                    scope.xdim === scope.ydim &&
                                    getDimType(scope.xdim) === 'ordinal') {

                                    return true;

                                } else {

                                    return false;
                                }
                            };

                            var drawNodes = function() {

                                if (isSameOrdDimGather()) {

                                    drawNodesForSameOrdDimGather();

                                } else {

                                    drawNodesForDifferentDim();
                                }


                            };

                            var drawNodesForDifferentDim = function() {

                                prepareScale();

                                calculateParametersOfNodes();

                                drawNodesInSVG();

                            }

                            var calculateParametersOfNodes = function() {

                                calculatePositionOfNodes();
                                calculateOffsetOfNodes();

                            };

                            var calculateParametersOfNodesForSameOrdDimGather = function() {

                                calculatePositionOfNodesForSameOrdDimGather();
                                calculateOffsetOfNodesForSameOrdDimGather();

                            };

                            var getKeyValue = function(dim) {

                                if (!dim) {
                                    return [''];
                                }

                                return dimSetting[dim].keyValue;
                            };

                            var getCalculatedPositions = function(dim) {

                                var keyValue = getKeyValue(dim);

                                var calculatedPosition = d3.map(keyValue)
                                    .entries()
                                    .map(function(d) {
                                        return +d.value.calculatedPosition;
                                    });

                                if (isDimTypeNumerical(getDimType(dim))) {

                                    calculatedPosition.sort(function(a, b) {
                                        return a - b;
                                    });

                                }

                                return calculatedPosition;


                            };



                            var getSortedIDs = function(dim) {

                                var keyValue = getKeyValue(dim);

                                var calculatedPosition = d3.map(keyValue)
                                    .entries()
                                    .map(function(d) {
                                        return +d.value.sortedID;
                                    });

                                return calculatedPosition;

                            };

                            //Returns Extents of dimension 
                            //              Scatter         Jitter      Gather
                            // ordinal      orig            orig        calculatedPoints
                            // semiordinal  SortedID        SortedID    calculatedPoints
                            // nominal      calculatedP     calculatedP calculatedPoints
                            var getExtent = function(dim) {

                                if (!dim) {

                                    return [-0.5, 0.5];
                                } else if (dim === 'Null') {

                                    return [-0.5, 0.5];

                                }

                                if (scope.config.isGather === 'gather') {

                                    if (dimSetting[dim].dimType === 'ordinal') {

                                        return getExtentFromOriginalExtent(dim);

                                    } else {

                                        return getExtentFromCalculatedPoints(dim);

                                    }
                                } else if (dimSetting[dim].dimType === 'ordinal') {

                                    return getExtentFromOriginalExtent(dim);

                                } else if (dimSetting[dim].dimType === 'semiOrdinal') {

                                    return getExtentFromSortedID(dim);

                                } else {

                                    return getExtentFromSortedID(dim);
                                }

                            };

                            var getDimType = function(dim) {

                                if (!dim) {
                                    return 'nominal';
                                } else if (dim === 'Null') {

                                    return 'nominal';


                                } else {

                                    return dimSetting[dim].dimType;
                                }
                            };

                            var getExtentFromSortedID = function(dim) {

                                var sortedID = getSortedIDs(dim);

                                var extent = d3.extent(sortedID);

                                return [extent[0] - marginForBorderOfAxis, extent[1] + marginForBorderOfAxis];

                            };


                            var getExtentFromCalculatedPoints = function(dim) {

                                calculatePositionOfCluster(dim);

                                var calculatedPoints = getCalculatedPositions(dim);

                                var max = calculatedPoints[calculatedPoints.length - 1];



                                var maxPadding = getLastIncrement(dim);

                                max = max + maxPadding;

                                return [0, max];




                            };

                            var getExtentFromCalculatedPointsForBinnedGather = function(dim) {

                                calculatePositionOfClusterForBinnedGather(dim);

                                var calculatedPoints = getCalculatedPositions(dim);

                                var max = calculatedPoints[calculatedPoints.length - 1];

                                return [0 - 0.5, max + 0.5];




                            };

                            var getLastIncrement = function(dim) {

                                if (!dim) {

                                    return;
                                }

                                var keyValue = dimSetting[dim].keyValue;
                                var increment;
                                var keyLength = d3.map(dimSetting[dim].keyValue).values().length;

                                var key = getKeyFromIndex(dim, keyLength - 1);

                                if (keyValue[key].isMinimized === true) {

                                    increment = marginClusterRatio;

                                } else {

                                    increment = 0.5;

                                }

                                return increment;

                            };

                            var getExtentFromOriginalExtent = function(dim) {

                                var originalValues = scope.data.map(function(d) {
                                    return +d[dim];
                                });

                                var extent = d3.extent(originalValues);

                                return [extent[0] - marginForBorderOfAxis, extent[1] + marginForBorderOfAxis];
                            };

                            var getExtentConsideringXY = function(xdim, ydim) {

                                var range = {};

                                var typeOfXYDim = findTypeOfXYDim();

                                var xRange, yRange;

                                if (typeOfXYDim === 'OrdOrd' && scope.config.isGather === 'gather') {

                                    doBinningAndSetKeys(xdim, scope.config.binSize);
                                    doBinningAndSetKeys(ydim, scope.config.binSize);

                                    xRange = getExtentFromCalculatedPoints(xdim);
                                    yRange = getExtentFromCalculatedPoints(ydim);


                                } else {

                                    xRange = getExtent(xdim);
                                    yRange = getExtent(ydim);

                                }
                                range.xRange = xRange;
                                range.yRange = yRange;

                                return range;

                            };

                            var prepareScale = function() {

                                var range = getExtentConsideringXY(scope.xdim, scope.ydim);

                                xRange = range.xRange;
                                yRange = range.yRange;



                                xScale = d3.scale.linear().range([0, width]);
                                xScale.domain(xRange);

                                xMap = function(d) {
                                    return xScale(xValue(d));
                                };

                                yScale = d3.scale.linear().range([height, 0]);
                                yScale.domain(yRange);
                                yMap = function(d) {
                                    return yScale(yValue(d));
                                };


                            };

                            function keyflip() {
                                shiftKey = d3.event.shiftKey || d3.event.metaKey;
                            }

                            function cross(a, b) {
                                var c = [],
                                    n = a.length,
                                    m = b.length,
                                    i, j;
                                for (i = -1; ++i < n;)
                                    for (j = -1; ++j < m;) c.push({
                                        x: a[i],
                                        i: i,
                                        y: b[j],
                                        j: j
                                    });
                                return c;
                            }

                            var restoreXYScaleForSameOrdDimGather = function() {

                                var xRange = getExtent(scope.xdim);
                                var yRange = getExtent(scope.ydim);

                                xScale = d3.scale.linear().range([0, width]);
                                xScale.domain(xRange);

                                xMap = function(d) {
                                    return xScale(xValue(d));
                                };

                                yScale = d3.scale.linear().range([height, 0]);
                                yScale.domain(yRange);
                                yMap = function(d) {
                                    return yScale(yValue(d));
                                };

                            };

                            var prepareScaleForSameOrdDimGather = function() {

                                var longAxisLength, shortAxisLength;

                                if (height < width) {

                                    longAxisLength = width;
                                    shortAxisLength = height;
                                } else {

                                    longAxisLength = height;
                                    shortAxisLength = width;
                                }

                                var virtualAxisLength = Math.sqrt(Math.pow(longAxisLength, 2) + Math.pow(shortAxisLength, 2));



                                var xRange = [0, 1];
                                var yRange = getExtent(scope.ydim);

                                xScale = d3.scale.linear().range([0, shortAxisLength]);
                                xScale.domain(xRange);

                                xMap = function(d) {
                                    return xScale(xValue(d));
                                };

                                yScale = d3.scale.linear().range([height, height - virtualAxisLength]);
                                yScale.domain(yRange);
                                yMap = function(d) {
                                    return yScale(yValue(d));
                                };

                            };

                            xOriginalValue = function(d) {

                                return d[scope.xdim];

                            };


                            yOriginalValue = function(d) {

                                return d[scope.ydim];
                            };



                            var dimOriginalValueConsideringBinning = function(dim) {

                                if (!dim) {

                                    return function(d) {
                                        return '';
                                    };
                                }

                                if (dimSetting[dim].isBinned) {

                                    return function(d) {

                                        return +dimSetting[dim].binnedData[d.id];
                                    };


                                } else {
                                    return function(d) {

                                        return d[dim];
                                    };
                                }
                            };



                            var colorOriginalValue = function(d) {

                                return d[scope.config.colorDim];
                            };



                            var calculatePositionOfNodes = function() {
                                //debugger;

                                if (scope.config.isGather === 'gather') {

                                    calculatePositionOfNodesForGather();

                                }

                                xValue = getPositionValueFunc(scope.xdim);
                                yValue = getPositionValueFunc(scope.ydim);


                            };

                            var calculatePositionOfNodesForSameOrdDimGather = function() {
                                //debugger;

                                var clusterSize = getClusterBox();
                                var range, height;


                                range = yScale.range();
                                height = range[0] - range[1];
                                getOptimalBinSize(scope.ydim, '', clusterSize.widthOfBox, height);

                                updateYScaleForSameOrdDimGather();
                                // calculatePositionOfCluster(scope.xdim);

                                xValue = getPositionValueFunc('');
                                yValue = getPositionValueFunc(scope.ydim);


                            };

                            var calculatePositionOfNodesForGather = function() {

                                var typeOfXYDim = findTypeOfXYDim();

                                if (typeOfXYDim === 'NomNom') {

                                    calculatePositionOfNodesForNomNomGather();

                                } else if (typeOfXYDim === 'OrdOrd') {

                                    calculatePositionOfNodesForOrdOrdGather();

                                } else {
                                    //Only one of them are ordinal -> binned gatherplot 

                                    calculatePositionOfNodesForBinnedGather();

                                }
                            };

                            var calculatePositionOfNodesForOrdOrdGather = function() {

                                var typeOfXYDim = findTypeOfXYDim();
                                var clusterSize = getClusterBox();
                                var range, height;

                                range = xScale.range();

                                calculatePositionOfCluster(scope.xdim);
                                calculatePositionOfCluster(scope.ydim);

                            };

                            var calculatePositionOfNodesForBinnedGather = function() {

                                var typeOfXYDim = findTypeOfXYDim();
                                var clusterSize = getClusterBox();
                                var range, height;

                                if (typeOfXYDim === 'XNomYOrd') {
                                    range = yScale.range();
                                    height = range[0] - range[1];
                                    getOptimalBinSize(scope.ydim, scope.xdim, clusterSize.widthOfBox, height);

                                    updateYScale();
                                    calculatePositionOfCluster(scope.xdim);
                                } else if (typeOfXYDim === 'XOrdYNom') {
                                    range = xScale.range();
                                    height = range[1] - range[0];
                                    getOptimalBinSize(scope.xdim, scope.ydim, clusterSize.heightOfBox, height);

                                    updateXScale();
                                    calculatePositionOfCluster(scope.ydim);
                                } else {


                                    calculatePositionOfCluster(scope.xdim);

                                    calculatePositionOfCluster(scope.ydim);


                                }

                            };



                            var updateYScale = function() {

                                var yRange = getExtentFromCalculatedPointsForBinnedGather(scope.ydim);

                                yScale = d3.scale.linear().range([height, 0]);
                                yScale.domain(yRange);
                                yMap = function(d) {
                                    return yScale(yValue(d));
                                };

                            };

                            var updateYScaleForSameOrdDimGather = function() {

                                var yRange = getExtentFromCalculatedPointsForBinnedGather(scope.ydim);

                                yScale.domain(yRange);
                                yMap = function(d) {
                                    return yScale(yValue(d));
                                };

                            };

                            var updateXScale = function() {

                                var xRange = getExtentFromCalculatedPointsForBinnedGather(scope.xdim);

                                xScale = d3.scale.linear().range([0, width]);
                                xScale.domain(xRange);
                                xMap = function(d) {
                                    return xScale(xValue(d));
                                };

                            };

                            var getOptimalBinSize = function(ordDim, nomDim, norDimLength, ordDimLength) {

                                var numBin = Math.floor(ordDimLength / maxDotSize);

                                var dotSize = maxDotSize;

                                var maxCrowdedBinCount = getMaxCrowdedBinCount(ordDim, nomDim, numBin);

                                var loopCount = 0;

                                var increment = numBin;
                                var previousIncrement = 1;

                                while (true) {


                                    if ((maxCrowdedBinCount) * dotSize > norDimLength) {

                                        increment = previousIncrement * 2;

                                    } else {

                                        increment = Math.round(previousIncrement * (-0.5));

                                    }

                                    numBin = numBin + increment;

                                    previousIncrement = Math.abs(increment);


                                    if (Math.abs(increment) < 2) {

                                        break;
                                    }

                                    maxCrowdedBinCount = getMaxCrowdedBinCount(ordDim, nomDim, numBin);
                                    dotSize = ordDimLength / numBin;

                                    loopCount = loopCount + 1;
                                }



                                console.log(loopCount + ": NumBin = " + numBin);

                                console.log(loopCount + ": increment = " + increment);


                                console.log(loopCount + ": maxCrowdedBinCount = " + getMaxCrowdedBinCount(ordDim, nomDim, numBin));

                                numBin = numBin + 1;

                                doBinningAndSetKeys(ordDim, numBin);

                            };

                            var getMaxCrowdedBinCount = function(ordDim, nomDim, binCount) {

                                var values = scope.data.map(function(d) {
                                    return +d[ordDim];
                                });

                                var ordinalScaleForGather = d3.scale.linear().domain(d3.extent(values));


                                var nestedData = d3.nest()
                                    .key(function(d) {
                                        return d[nomDim];
                                    })
                                    .entries(scope.data);

                                var maxValues = nestedData.map(function(d) {

                                    var values = d.values.map(function(d) {
                                        return +d[ordDim];
                                    });

                                    var data = d3.layout.histogram()
                                        .bins(binCount)
                                        (values);

                                    // console.log(data.bins());

                                    return d3.max(data, function(d) {
                                        return +d.y;
                                    });
                                });

                                return d3.max(maxValues) + 1;

                            }

                            var findTypeOfXYDim = function() {

                                if (scope.xdim === 'Null') {

                                    scope.xdim = '';
                                    scope.config.xdim = '';
                                }

                                if (scope.ydim === 'Null') {

                                    scope.ydim = '';
                                    scope.config.ydim = '';
                                }

                                if (scope.config.colorDim === 'Null') {

                                    scope.config.colorDim = '';
                                }


                                var xDimType = getDimType(scope.xdim);
                                var yDimType = getDimType(scope.ydim);

                                if (xDimType === 'ordinal' && yDimType === 'ordinal') {

                                    if (scope.xdim === scope.ydim) {
                                        return 'SameOrd';
                                    } else {
                                        return 'OrdOrd';
                                    }
                                } else if (xDimType !== 'ordinal' && yDimType !== 'ordinal') {
                                    return 'NomNom';
                                } else if (xDimType === 'ordinal' && yDimType !== 'ordinal') {
                                    return 'XOrdYNom';
                                } else if (xDimType !== 'ordinal' && yDimType === 'ordinal') {
                                    return 'XNomYOrd';
                                }

                            };

                            var calculatePositionOfNodesForNomNomGather = function() {

                                calculatePositionOfCluster(scope.xdim);
                                calculatePositionOfCluster(scope.ydim);

                            }

                            var calculatePositionOfCluster = function(dim) {

                                if (!dim) {

                                    return;
                                } else if (dim === 'Null') {
                                    return;
                                }

                                var keyValue = dimSetting[dim].keyValue;
                                var increment;
                                var previousIncrement;



                                d3.map(keyValue).entries().forEach(function(d, i, all) {


                                    if (i === 0) {
                                        if (d.value.isMinimized === true) {

                                            d.value.calculatedPosition = marginClusterRatio;

                                        } else {

                                            d.value.calculatedPosition = 0.5;

                                        }

                                        d.value.calculatedPosition = d.value.calculatedPosition;

                                        return;
                                    }

                                    if (all[i - 1].value.isMinimized === true) {

                                        previousIncrement = marginClusterRatio;

                                    } else {

                                        previousIncrement = 0.5;

                                    }


                                    if (d.value.isMinimized === true) {

                                        increment = marginClusterRatio;

                                    } else {

                                        increment = 0.5;

                                    }

                                    d.value.calculatedPosition = all[i - 1].value.calculatedPosition + increment + previousIncrement;

                                });

                            };

                            var calculatePositionOfClusterForBinnedGather = function(dim) {


                                var keyValue = dimSetting[dim].keyValue;

                                var keys = Object.keys(keyValue);

                                keys.sort(function(a, b) {
                                    return a - b;
                                });

                                for (var i = 0; i < keys.length; i++) {

                                    keyValue[keys[i]].value = {};

                                    keyValue[keys[i]].value.calculatedPosition = i + 0.5;
                                }

                            };


                            var getPositionValueFunc = function(dimName) {

                                if (!dimName) {

                                    return function(d) {
                                        return 0;
                                    };
                                }

                                var dimType = dimSetting[dimName].dimType;
                                var dimNameClosure = dimName;

                                // Follow the dimValue branch logic
                                //                  Scatter         Jitter       Gather
                                // nominal           sortedID      sortedID       calculatedID
                                // SemiOrdi         sortedID       sortedID       calculatedID
                                // ordinal           orig          orig           calculatedIDFromBin

                                var calculatedPositionValueFunc = function(d) {
                                    return dimSetting[dimNameClosure].keyValue[d[dimNameClosure]].calculatedPosition;
                                };

                                var origValueFunc = function(d) {

                                    return +d[dimNameClosure];
                                };

                                var calculatedPositionWithBinValueFunc = function(d) {
                                    var binKey = +dimSetting[dimNameClosure].binnedData[d.id];
                                    if (!dimSetting[dimNameClosure].keyValue[binKey]) {

                                        //console.log(binKey);
                                    }

                                    var positionWithBinKey = dimSetting[dimNameClosure].keyValue[binKey].calculatedPosition;

                                    return +positionWithBinKey;
                                };

                                var sortedPositionValueFunc = function(d) {
                                    return dimSetting[dimNameClosure].keyValue[d[dimNameClosure]].sortedID;
                                };

                                if (dimType === 'nominal') {

                                    if (scope.config.isGather === 'gather') {

                                        return calculatedPositionValueFunc;
                                    } else {

                                        return sortedPositionValueFunc;
                                    }

                                } else if (dimType === 'semiOrdinal') {

                                    if (scope.config.isGather === 'gather') {

                                        return calculatedPositionValueFunc;

                                    } else {

                                        return sortedPositionValueFunc;
                                    }

                                } else if (dimType === 'ordinal') {


                                    if (scope.config.isGather === 'gather') {
                                        return calculatedPositionWithBinValueFunc;

                                    } else {
                                        return origValueFunc;
                                    }

                                } else {

                                    console.log("Unsupported DimName in getDimValueFunc");
                                }

                            };



                            var calculateOffsetOfNodes = function() {

                                if (scope.config.isGather === 'scatter') {

                                    setOffsetOfNodesForScatter();

                                } else if (scope.config.isGather === 'jitter') {

                                    setOffsetOfNodesForJitter();

                                } else if (scope.config.isGather === 'gather') {


                                    setOffsetOfNodesForGather();

                                }

                            };

                            var calculateOffsetOfNodesForSameOrdDimGather = function() {


                                setOffsetOfNodesForGatherForSameOrdDimGather();


                            };

                            var setOffsetOfNodesForScatter = function() {

                                scope.data.forEach(function(d) {

                                    d.XOffset = 0;
                                    d.YOffset = 0;

                                });

                                assignSizeOfNodesForScatterAndJitter();

                            };

                            var assignSizeOfNodesForScatterAndJitter = function() {



                                scope.data.forEach(function(d) {

                                    d.nodeWidth = maxDotSize;
                                    d.nodeHeight = maxDotSize;

                                });
                            };

                            var setOffsetOfNodesForJitter = function() {


                                var SDforJitter = getSDforJitter();

                                var xNormalGenerator = d3.random.normal(0, SDforJitter.xSD);
                                var yNormalGenerator = d3.random.normal(0, SDforJitter.ySD);

                                scope.data.forEach(function(d) {


                                    d.XOffset = xNormalGenerator();
                                    d.YOffset = yNormalGenerator();

                                });

                                assignSizeOfNodesForScatterAndJitter();

                            };

                            var setOffsetOfNodesForGather = function() {

                                makeNestedData();

                                assignClusterIDOfNodes();
                                updateClusterSizeInNestedData();
                                getNodesSizeAndOffsetPosition();
                                // assignOffsetForGather();

                            };

                            var setOffsetOfNodesForGatherForSameOrdDimGather = function() {

                                makeNestedDataForSameOrdDimGather();

                                assignClusterIDOfNodes();
                                updateClusterSizeInNestedData();
                                getNodesSizeAndOffsetPosition();
                                // assignOffsetForGather();

                            };

                            var makeNestedData = function() {


                                // debugger;

                                var xOriginalValueWithBinning = dimOriginalValueConsideringBinning(scope.xdim);

                                var yOriginalValueWithBinning = dimOriginalValueConsideringBinning(scope.ydim);

                                nest = d3.nest()
                                    .key(xOriginalValueWithBinning)
                                    .key(yOriginalValueWithBinning)
                                    .sortValues(sortFuncByColorDimension())
                                    .entries(scope.data);


                            };

                            var makeNestedDataForSameOrdDimGather = function() {


                                // debugger;

                                var xOriginalValueWithBinning = dimOriginalValueConsideringBinning('');

                                var yOriginalValueWithBinning = dimOriginalValueConsideringBinning(scope.ydim);

                                nest = d3.nest()
                                    .key(xOriginalValueWithBinning)
                                    .key(yOriginalValueWithBinning)
                                    .sortValues(sortFuncByColorDimension())
                                    .entries(scope.data);


                            };

                            var assignClusterIDOfNodes = function() {

                                assignClusterIDOfNodesInTwoKeyNestedItems(nest);

                            };

                            var assignClusterIDOfNodesInTwoKeyNestedItems = function(nest) {

                                nest.forEach(function(d, i, j) {

                                    d.values.forEach(function(d, i, j) {

                                        d.values.forEach(function(d, i, j) {

                                            d.clusterID = i;

                                        });

                                    });

                                });


                            };

                            var assignClusterIDOfNodesInOneKeyNestedItems = function(nest) {

                                nest.forEach(function(d, i, j) {

                                    d.values.forEach(function(d, i, j) {

                                        d.clusterID = i;


                                    });

                                });


                            };

                            var updateClusterSizeInNestedData = function() {

                                nest.forEach(function(d, i, j) {

                                    d.values.forEach(function(d, i, j) {

                                        d.numOfElement = d.values.length;

                                    });

                                });


                            };

                            var sortFuncByColorDimension = function() {

                                var colorDim = scope.config.colorDim;

                                if (!colorDim) {
                                    return function(a, b) {
                                        return a;
                                    };
                                } else {

                                    // debugger;

                                    if (isDimTypeNumerical(dimSetting[colorDim].dimType)) {

                                        return numericalDimSortFunc(colorDim);

                                    } else {

                                        return nominalDimSortFunc(colorDim);

                                    }


                                }

                            };

                            var nominalDimSortFunc = function(dim) {

                                var tempDimSetting = dimSetting[dim];

                                return function(a, b) {
                                    var myDim = dim;
                                    return tempDimSetting.keyValue[a[myDim]].sortedID - tempDimSetting.keyValue[b[myDim]].sortedID;
                                };

                            };

                            var numericalDimSortFunc = function(dim) {

                                return function(a, b) {
                                    return a[dim] - b[dim];
                                };
                            };

                            var isDimTypeNumerical = function(dimType) {

                                if (dimType === 'nominal') {

                                    return false;

                                } else if (dimType === 'ordinal' || dimType === 'semiOrdinal') {

                                    return true;
                                } else {

                                    alert("Unidentified dimension type");
                                }
                            };

                            var getClusterBox = function() {

                                var Xmargin, Ymargin;
                                var typeOfXYDim = findTypeOfXYDim();

                                if (typeOfXYDim === 'NomNom') {

                                    Xmargin = marginClusterRatio;
                                    Ymargin = marginClusterRatio;
                                } else if (typeOfXYDim === 'XNomYOrd') {

                                    Xmargin = marginClusterRatio;
                                    Ymargin = 0;
                                } else if (typeOfXYDim === 'XOrdYNom') {

                                    Xmargin = 0;
                                    Ymargin = marginClusterRatio;
                                } else if (typeOfXYDim === 'OrdOrd') {

                                    Xmargin = marginClusterRatio;
                                    Ymargin = marginClusterRatio;

                                } else {

                                    Xmargin = 0;
                                    Ymargin = 0;
                                }


                                return {
                                    widthOfBox: xScale(1 - 2 * Xmargin) - xScale(0),
                                    heightOfBox: yScale(0) - yScale(1 - 2 * Ymargin)
                                };

                            };


                            var getNodesSizeForAbsolute = function() {

                                var maxNumElementInCluster = getClusterWithMaximumPopulation();
                                // console.log('maxNumElementInCluster = ' + maxNumElementInCluster )
                                var box = getClusterBox();
                                var size = calculateNodesSizeForAbsolute(box, maxNumElementInCluster);

                                return size;

                            };


                            var getNodesSizeAndOffsetPosition = function() {

                                nest.forEach(function(d, i, j) {

                                    var xKey = d.key;

                                    d.values.forEach(function(d, i, j) {

                                        var yKey = d.key;

                                        assignNodesOffsetByCluster(d.values, xKey, yKey);

                                    });

                                });


                            };


                            var assignNodesOffsetByCluster = function(cluster, xKey, yKey) {

                                var box = getClusterBox();

                                assignNodesOffsetConsideringAspectRatio(cluster, box)


                                updateNodesOffsetForMinimized(cluster, xKey, yKey);
                                updateNodesSizeForMinimized(cluster, xKey, yKey);

                            };

                            var assignNodesOffsetConsideringAspectRatio = function(cluster, box) {

                                if (box.widthOfBox > box.heightOfBox) {

                                    assignNodesOffsetHorizontallyByCluster(cluster, box);

                                } else {

                                    assignNodesOffsetVerticallyByCluster(cluster, box);
                                }



                            };

                            var updateNodesSizeForMinimized = function(cluster, xKey, yKey) {

                                if (isMinimized(scope.xdim, xKey)) {

                                    makeAbsoluteSize(cluster, 'nodeWidth');
                                }

                                if (isMinimized(scope.ydim, yKey)) {

                                    makeAbsoluteSize(cluster, 'nodeHeight');
                                }

                            };

                            var updateNodesOffsetForMinimized = function(cluster, xKey, yKey) {

                                if (isMinimized(scope.xdim, xKey)) {

                                    makeZeroOffset(cluster, 'XOffset');

                                }

                                if (isMinimized(scope.ydim, yKey)) {

                                    makeZeroOffset(cluster, 'YOffset');
                                }



                            };

                            var isMinimized = function(dim, key) {

                                if (!dim) {

                                    return false;
                                }

                                if (!key) {

                                    return false;
                                }

                                if (!scope.config.isInteractiveAxis) {

                                    return false;
                                }

                                return (dimSetting[dim].keyValue[key].isMinimized);
                            };

                            var makeZeroOffset = function(cluster, offset) {

                                cluster.forEach(function(d) {

                                    d[offset] = 0;

                                });
                            };
                            var makeAbsoluteSize = function(cluster, nodeSize) {

                                var absoulteSize = getNodesSizeForAbsolute();

                                cluster.forEach(function(d) {

                                    d[nodeSize] = absoulteSize;

                                });
                            };




                            var assignNodesOffsetLongShortEdge = function(longEdge, shortEdge, cluster) {

                                var numElement = getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping(longEdge, shortEdge, cluster.length);
                                if (isThemeRiverCondition(longEdge, shortEdge, numElement)) {

                                    numElement = getNumOfElementForThemeRiver(longEdge, shortEdge, cluster.length);
                                    if (numElement.numElementInShortEdge === 2) {
                                        console.log('Hey');
                                    }
                                }
                                var nodeSize = getNodeSizeAbsoluteOrRelative(longEdge, shortEdge, numElement.numElementInLongEdge, numElement.numElementInShortEdge);
                                var offsetForCenterPosition = calculateOffsetForCenterPosition(nodeSize.lengthInLongEdge, nodeSize.lengthInShortEdge, numElement.numElementInLongEdge, numElement.numElementInShortEdge);


                                return {
                                    numElement: numElement,
                                    nodeSize: nodeSize,
                                    offsetForCenterPosition: offsetForCenterPosition
                                };


                            };

                            var assignNodesOffsetLongShortEdgeLens = function(longEdge, shortEdge, cluster, size) {

                                var numElement = getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping(longEdge, shortEdge, cluster.length, maxNum);
                                if (isThemeRiverCondition(longEdge, shortEdge, numElement)) {

                                    numElement = getNumOfElementForThemeRiver(longEdge, shortEdge, cluster.length);
                                }
                                var nodeSize = getNodeSizeAbsoluteOrRelative(longEdge, shortEdge, numElement.numElementInLongEdge, numElement.numElementInShortEdge);
                                var offsetForCenterPosition = calculateOffsetForCenterPosition(nodeSize.lengthInLongEdge, nodeSize.lengthInShortEdge, numElement.numElementInLongEdge, numElement.numElementInShortEdge);


                                return {
                                    numElement: numElement,
                                    nodeSize: nodeSize,
                                    offsetForCenterPosition: offsetForCenterPosition
                                };


                            };

                            var isThemeRiverCondition = function(longEdge, shortEdge, numElement) {

                                if (longEdge / shortEdge > 3) {

                                    return true;
                                } else {

                                    return false;
                                }
                            };

                            var getNumOfElementForThemeRiver = function(longEdge, shortEdge, numElement) {

                                var numElementInShortEdge = Math.ceil(shortEdge / getNodesSizeForAbsolute());


                                if (numElementInShortEdge == 2 && getNodesSizeForAbsolute() < 1) {

                                    numElementInShortEdge = 1;


                                }

                                var numElementInLongEdge = Math.ceil(numElement / numElementInShortEdge);


                                return {
                                    numElementInShortEdge: numElementInShortEdge,
                                    numElementInLongEdge: numElementInLongEdge
                                };


                            };

                            var getNodeSizeAbsoluteOrRelative = function(longEdge, shortEdge, numElementInLongEdge, numElementInShortEdge) {

                                var lengthInLongEdge, lengthInShortEdge;

                                if (scope.config.relativeMode === "absolute") {

                                    lengthInLongEdge = getNodesSizeForAbsolute();
                                    lengthInShortEdge = lengthInLongEdge;

                                } else {
                                    lengthInLongEdge = longEdge / numElementInLongEdge;
                                    lengthInShortEdge = shortEdge / numElementInShortEdge;
                                }

                                return {
                                    lengthInLongEdge: lengthInLongEdge,
                                    lengthInShortEdge: lengthInShortEdge
                                };

                            };

                            var handleOffsetRectLensHorizontally = function(cluster, box, size) {

                                var nodeHeight = size;
                                var nodeWidth = size;

                                var numOfElement = getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping(box.widthOfBox, box.heightOfBox, cluster.length);
                                var numElementInShortEdge = numOfElement.numElementInShortEdge;
                                var numElementInLongEdge = numOfElement.numElementInLongEdge;
                                var offsetInShortEdge = nodeHeight * numElementInShortEdge / 2;
                                var offsetInLongEdge = nodeWidth * numElementInLongEdge / 2;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeWidthLens = nodeWidth;
                                    d.nodeHeightLens = nodeHeight;




                                    d.YOffsetLens = (d.clusterID % numElementInShortEdge) * nodeHeight - offsetInShortEdge + nodeHeight;
                                    d.XOffsetLens = Math.floor(d.clusterID / numElementInShortEdge) * nodeWidth - offsetInLongEdge;

                                });

                            };

                            var handleOffsetRectLensVertically = function(cluster, box, size) {

                                var nodeHeight = size;
                                var nodeWidth = size;

                                var numOfElement = getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping(box.heightOfBox, box.widthOfBox, cluster.length);
                                var numElementInShortEdge = numOfElement.numElementInShortEdge;
                                var numElementInLongEdge = numOfElement.numElementInLongEdge;
                                var offsetInShortEdge = nodeWidth * numElementInShortEdge / 2;
                                var offsetInLongEdge = nodeHeight * numElementInLongEdge / 2;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeHeightLens = nodeHeight;
                                    d.nodeWidthLens = nodeWidth;

                                    d.XOffsetLens = (d.clusterID % numElementInShortEdge) * nodeWidth - offsetInShortEdge;
                                    d.YOffsetLens = Math.floor(d.clusterID / numElementInShortEdge) * nodeHeight - offsetInLongEdge + nodeHeight;

                                });

                            };

                            var handleOffsetHistLensHorizontally = function(cluster, box, size) {

                                var nodeHeight = size;
                                var nodeWidth = size;
                                var numElementInShortEdge = Math.round(box.heightOfBox / size);
                                var numElementInLongEdge = Math.round(box.widthOfBox / size);
                                var offsetInShortEdge = nodeHeight * numElementInShortEdge / 2;
                                var offsetInLongEdge = nodeWidth * numElementInLongEdge / 2;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeWidthLens = nodeWidth;
                                    d.nodeHeightLens = nodeHeight;

                                    d.YOffsetLens = (d.clusterID % numElementInShortEdge) * nodeHeight - offsetInShortEdge + nodeHeight;
                                    d.XOffsetLens = Math.floor(d.clusterID / numElementInShortEdge) * nodeWidth - offsetInLongEdge;

                                });

                            };

                            var handleOffsetHistLensVertically = function(cluster, box, size) {

                                var nodeHeight = size;
                                var nodeWidth = size;
                                var numElementInShortEdge = Math.round(box.widthOfBox / size);
                                var numElementInLongEdge = Math.round(box.heightOfBox / size);
                                var offsetInShortEdge = nodeHeight * numElementInShortEdge / 2;
                                var offsetInLongEdge = nodeWidth * numElementInLongEdge / 2;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeWidthLens = nodeWidth;
                                    d.nodeHeightLens = nodeHeight;

                                    d.XOffsetLens = (d.clusterID % numElementInShortEdge) * nodeWidth - offsetInShortEdge;
                                    d.YOffsetLens = Math.floor(d.clusterID / numElementInShortEdge) * nodeHeight - offsetInLongEdge + nodeHeight;

                                });

                            };


                            var assignNodesOffsetHorizontallyByCluster = function(cluster, box) {

                                var offsetAndSizeInfo = assignNodesOffsetLongShortEdge(box.widthOfBox, box.heightOfBox, cluster);

                                var nodeHeight = offsetAndSizeInfo.nodeSize.lengthInShortEdge;
                                var nodeWidth = offsetAndSizeInfo.nodeSize.lengthInLongEdge;
                                var numElementInShortEdge = offsetAndSizeInfo.numElement.numElementInShortEdge;
                                var numElementInLongEdge = offsetAndSizeInfo.numElement.numElementInLongEdge;
                                var offsetInShortEdge = offsetAndSizeInfo.offsetForCenterPosition.offsetInShortEdge;
                                var offsetInLongEdge = offsetAndSizeInfo.offsetForCenterPosition.offsetInLongEdge;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeWidth = nodeWidth;
                                    d.nodeHeight = nodeHeight;




                                    d.YOffset = (d.clusterID % numElementInShortEdge) * nodeHeight - offsetInShortEdge + nodeHeight;
                                    d.XOffset = Math.floor(d.clusterID / numElementInShortEdge) * nodeWidth - offsetInLongEdge;

                                });

                            };

                            var assignNodesOffsetVerticallyByCluster = function(cluster, box) {

                                var offsetAndSizeInfo = assignNodesOffsetLongShortEdge(box.heightOfBox, box.widthOfBox, cluster);

                                var nodeHeight = offsetAndSizeInfo.nodeSize.lengthInLongEdge;
                                var nodeWidth = offsetAndSizeInfo.nodeSize.lengthInShortEdge;
                                var numElementInShortEdge = offsetAndSizeInfo.numElement.numElementInShortEdge;
                                var numElementInLongEdge = offsetAndSizeInfo.numElement.numElementInLongEdge;
                                var offsetInShortEdge = offsetAndSizeInfo.offsetForCenterPosition.offsetInShortEdge;
                                var offsetInLongEdge = offsetAndSizeInfo.offsetForCenterPosition.offsetInLongEdge;

                                cluster.forEach(function(d, i, j) {

                                    d.nodeHeight = nodeHeight;
                                    d.nodeWidth = nodeWidth;

                                    d.XOffset = (d.clusterID % numElementInShortEdge) * nodeWidth - offsetInShortEdge;
                                    d.YOffset = Math.floor(d.clusterID / numElementInShortEdge) * nodeHeight - offsetInLongEdge + nodeHeight;

                                });

                            };

                            var calculateOffsetForCenterPosition = function(nodeLengthInLongEdge, nodeLengthInShortEdge, numElementInLongEdge, numElementInShortEdge) {

                                var offsetInShortEdgeForCenterPosition;
                                var offsetInLongEdgeForCenterPosition;

                                offsetInShortEdgeForCenterPosition = numElementInShortEdge * nodeLengthInShortEdge / 2;
                                offsetInLongEdgeForCenterPosition = numElementInLongEdge * nodeLengthInLongEdge / 2;

                                return {
                                    offsetInShortEdge: offsetInShortEdgeForCenterPosition,
                                    offsetInLongEdge: offsetInLongEdgeForCenterPosition
                                };
                            };

                            var getClusterWithMaximumPopulation = function() {

                                return getClusterWithMaximumPopulationFromTwoKeyNestedItems(nest);
                            };

                            var getClusterWithMaximumPopulationFromTwoKeyNestedItems = function(nest) {

                                return d3.max(nest, function(d) {

                                    return d3.max(d.values, function(d) {

                                        return d.numOfElement;
                                    });
                                });

                            };

                            var getClusterWithMaximumPopulationFromOneKeyNestedItems = function(nest) {

                                return d3.max(nest, function(d) {

                                    return d.values.length;
                                });

                            };

                            var calculateNodesSizeForAbsolute = function(box, maxNumber) {

                                if (box.widthOfBox > box.heightOfBox) {

                                    return calculateNodesSizeWithLongAndShortEdges(box.widthOfBox, box.heightOfBox, maxNumber);

                                } else {

                                    return calculateNodesSizeWithLongAndShortEdges(box.heightOfBox, box.widthOfBox, maxNumber);
                                }
                            };

                            var calculateNodesSizeWithLongAndShortEdges = function(longEdge, shortEdge, number) {


                                var numElement = getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping(longEdge, shortEdge, number);

                                return shortEdge / numElement.numElementInShortEdge;

                            };

                            var getNumOfElementInLongAndShortEdgeUsingAspectRatioKeeping = function(longEdge, shortEdge, number) {

                                var numElementInShortEdge = 0,
                                    numElementInLongEdge,
                                    sizeNode, lengthCandidate;



                                do {

                                    numElementInShortEdge++;
                                    sizeNode = shortEdge / numElementInShortEdge;
                                    lengthCandidate = sizeNode * number / numElementInShortEdge;

                                } while (lengthCandidate > longEdge);

                                numElementInLongEdge = Math.ceil(number / numElementInShortEdge);

                                return {
                                    numElementInShortEdge: numElementInShortEdge,
                                    numElementInLongEdge: numElementInLongEdge
                                };


                            };



                            var getSDforJitter = function() {

                                var nominalBox = getClusterBox();
                                var probFactor = 0.15;

                                var xSD = nominalBox.widthOfBox * probFactor;
                                var ySD = nominalBox.heightOfBox * probFactor;

                                return {
                                    xSD: xSD,
                                    ySD: ySD
                                };

                            };



                            var drawNodesInSVG = function() {

                                getColorOfNodes();
                                getShapeOfNodes();
                                writeNodesInSVG();


                            };

                            var drawNodesInSVGForSameOrdDimGather = function() {

                                getColorOfNodes();
                                getShapeOfNodes();
                                writeNodesInSVGForSameOrdDimGather();


                            };

                            var getColorOfNodes = function() {

                                if (!scope.config.colorDim) {
                                    color = colorNominal;
                                    return;
                                }

                                if (dimSetting[scope.config.colorDim].dimType === 'ordinal') {

                                    var colorDomain = d3.extent(scope.data, function(d) {
                                        return +d[scope.config.colorDim];
                                    });

                                    colorScaleForHeatMap = d3.scale.linear()
                                        .range(["#98c8fd", "08306b"])
                                        .domain(colorDomain)
                                        .interpolate(d3.interpolateHsl);

                                    color = colorScaleForHeatMap;
                                } else {

                                    color = colorNominal;
                                }

                            };

                            var getShapeOfNodes = function() {

                            };

                            var writeNodesInSVG = function() {
                                nodeGroup.selectAll("text").remove();
                                // debugger;

                                // nodeGroup.attr("transform", "translate(" + margin + "," + margin + ") rotate(0 80 660)");

                                nodeGroup.attr("transform", "translate(0,0) rotate(0 80 660)");

                                var clusterInfo = {};

                                nodeGroup.selectAll(".dot")
                                    .data(scope.data, function(d) {
                                        return +d.id;
                                    })
                                    .style("fill", function(d) {
                                        return color(d[scope.config.colorDim]-1);
                                    })
                                    .transition()
                                    .duration(1500)
                                    .attr("x", function(d) {
                                        var clusterIndex = d.cluster
                                        if(clusterInfo[clusterIndex]==undefined) {
                                            clusterInfo[clusterIndex] = {};
                                            clusterInfo[clusterIndex].index = clusterIndex;
                                            clusterInfo[clusterIndex].num = 1;
                                            clusterInfo[clusterIndex].X = parseFloat(d.X);
                                            clusterInfo[clusterIndex].Y = parseFloat(d.Y);
                                            clusterInfo[clusterIndex].keywords = [
                                                d['cluster top keyword 1'],
                                                d['cluster top keyword 2'],
                                                d['cluster top keyword 3'],
                                                d['cluster top keyword 4'],
                                                d['cluster top keyword 5']
                                            ];
                                        } else {
                                            clusterInfo[clusterIndex].num += 1;
                                            clusterInfo[clusterIndex].X += parseFloat(d.X);
                                            clusterInfo[clusterIndex].Y += parseFloat(d.Y);
                                        }
                                        return xMap(d);
                                    })
                                    .attr("y", yMap)
                                    .attr("width", function(d) {
                                        // console.log(initialSquareLenth);
                                        return +d.nodeWidth*1.5;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeight*1.5;
                                    })
                                    .attr("rx", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .attr("ry", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .attr("transform", function(d, i) {

                                        // if (d.cancer== "Cancer") {
                                        //     console.log(height);
                                        // }
                                        return "translate(" + (d.XOffset) + "," + (-(d.YOffset)) + ") ";
                                    });


                                clusterInfo = $.map(clusterInfo, function(d,i) { return [d]; });
                                for(var i=0;i<clusterInfo.length;i++) {
                                    var num = clusterInfo[i].num;
                                    clusterInfo[i].X/=num;
                                    clusterInfo[i].Y/=num;
                                }
                                
                                nodeGroup.selectAll(".text")
                                    .data(clusterInfo).enter()
                                    .append("text")
                                    .attr('x', xMap)
                                    .attr('y', yMap)
                                    .style("fill", "black")
                                    .text(function(d){
                                        return d.keywords.join(' ');
                                    })
                                    .attr('font-family','sans-serif')
                                    .attr('text-anchor','middle')
                                    .style('fill','black')
                                    .style('pointer-events', 'none');

                            };

                            var writeNodesInSVGForSameOrdDimGather = function() {
                                // debugger;

                                nodeGroup.selectAll(".dot")
                                    .data(scope.data, function(d) {
                                        return +d.id;
                                    })
                                    .style("fill", function(d) {
                                        return color(d[scope.config.colorDim]);
                                    })
                                    .transition()
                                    .duration(0)
                                    .attr("x", xMap)
                                    .attr("y", yMap)
                                    .attr("width", function(d) {
                                        // console.log(initialSquareLenth);
                                        return +d.nodeWidth;
                                    })
                                    .attr("height", function(d) {
                                        return +d.nodeHeight;
                                    })
                                    .attr("rx", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .attr("ry", function(d) {
                                        return scope.round ? +5 : 0;
                                    })
                                    .attr("transform", function(d, i) {

                                        // if (d.cancer== "Cancer") {
                                        //     console.log(height);
                                        // }
                                        return "translate(" + (d.XOffset) + "," + (-(d.YOffset)) + ") ";
                                    });

                                var angleRad = Math.atan(height / width);

                                var angleDeg = 90 - angleRad * 180 / Math.PI;


                                nodeGroup.attr("transform", " translate(" + margin + "," + margin + ")  rotate(" + angleDeg + "," + "0" + "," + yScale.range()[0] + ")");

                            };

                            var labelGenerator = function(dimName) {

                                if (!dimName) {

                                    return function(d) {
                                        return '';
                                    };
                                } else if ((dimSetting[dimName].dimType === 'ordinal')) {

                                    return function(d, i) {

                                        return +d;
                                    };
                                } else if ((dimSetting[dimName].dimType === 'semiOrdinal')) {

                                    return function(d, i) {

                                        return d3.map(dimSetting[dimName].keyValue).keys()[i];
                                    };
                                } else {

                                    return function(d) {



                                        return getKeys(dimName)[d];

                                    };
                                }


                            };

                            var labelGeneratorForGather = function(dimName) {

                                if (!dimName) {

                                    return function(d) {
                                        return '';
                                    };
                                } else if (dimSetting[dimName].dimType === 'ordinal') {

                                    var binDistanceFormatter = d3.format("3,.1f");

                                    return function(d, i) {

                                        var binValue = d3.map(dimSetting[dimName].keyValue).keys()[i];

                                        return binDistanceFormatter(+binValue) + '\u00B1' + binDistanceFormatter(+dimSetting[dimName].halfOfBinDistance);
                                    };
                                } else if (dimSetting[dimName].dimType === 'semiOrdinal') {

                                    return function(d, i) {

                                        return d3.map(dimSetting[dimName].keyValue).keys()[i];
                                    };
                                } else {

                                    return function(d, i) {



                                        return getKeys(dimName)[i];

                                    };
                                }


                            };

                            var labelGeneratorForOrdinalGather = function(dim) {

                                var keyValue = dimSetting[dim].keyValue;

                                var keys = Object.keys(keyValue)
                                    .sort(function(a, b) {
                                        return a - b;
                                    });

                                var binDistanceFormatter = d3.format("3,.0f");


                                return function(d, i) {

                                    return binDistanceFormatter(+keys[d]);

                                };


                            };

                            var tickGenerator = function(dimName) {

                                if (!dimName) {
                                    return 0;
                                } else if (dimSetting[dimName].dimType === 'ordinal') {

                                    return 8;

                                } else {

                                    return getKeys(dimName).length;
                                }
                            };

                            var tickValueGeneratorForGather = function(dimName) {

                                if (!dimName) {
                                    return [];

                                }
                                return getCalculatedPositions(dimName);

                            };

                            var tickValueGeneratorForSameOrdGather = function(dimName) {

                                if (!dimName) {
                                    return [];

                                }


                                var originalPositions = getCalculatedPositions(dimName);


                                var samplingRate = getSamplingRateForOrdinalGather(dimName);

                                var sampledPositions = originalPositions.filter(function(d, i) {
                                    return (i % samplingRate === 0);
                                });



                                sampledPositions = sampledPositions.map(function(d) {
                                    return d + Math.floor(samplingRate / 0.5);
                                })

                                sampledPositions.pop();

                                return sampledPositions;

                            };

                            var tickValueGeneratorForOrdinalGather = function(dimName) {

                                if (!dimName) {
                                    return [];

                                }


                                var originalPositions = getCalculatedPositions(dimName);


                                var samplingRate = getSamplingRateForOrdinalGather(dimName);

                                var sampledPositions = originalPositions.filter(function(d, i) {
                                    return (i % samplingRate === 0);
                                });



                                sampledPositions = sampledPositions.map(function(d) {
                                    return d + Math.floor(samplingRate / 2);
                                })

                                sampledPositions.pop();

                                return sampledPositions;

                            };

                            var getSamplingRateForOrdinalGather = function(dimName) {

                                var originalPositions = getCalculatedPositions(dimName);

                                var dimLength = originalPositions.length;

                                return Math.floor(dimLength / 7);

                            }


                            var drawAxes = function() {

                                if (isSameOrdDimGather()) {

                                    drawAxesForSameOrdDimGather();
                                } else {

                                    drawAxesForDifferentDim();
                                }

                            };

                            var drawAxesForDifferentDim = function() {

                                drawAxesLinesAndTicks();
                                drawAxesLabel();

                            }

                            var drawAxesForSameOrdDimGather = function() {

                                restoreXYScaleForSameOrdDimGather();

                                drawAxesLinesAndTicksForSameOrdDimGather();
                                drawAxesLabel();
                                setStylesForAxesAndTicks();

                            };

                            var drawAxesLinesAndTicks = function() {

                                if (scope.config.isGather === 'gather') {

                                    drawAxesLinesAndTicksForGather();

                                } else {

                                    drawAxesLinesAndTicksForScatter();
                                }

                                setStylesForAxesAndTicks();


                            };

                            var setStylesForAxesAndTicks = function() {

                                svg.selectAll(".domain")
                                    .style("stroke", "black")
                                    .style("stroke-width", 1)
                                    .style("fill", "none");

                                svg.selectAll(".bracket")
                                    .style("stroke", "black")
                                    .style("stroke-width", 1)
                                    .style("fill", "none");


                            };

                            var drawAxesLinesAndTicksForScatter = function() {

                                svg.selectAll(".axis").remove();

                                drawXAxisLinesAndTicksForScatter();
                                drawYAxisLinesAndTicksForScatter();

                            };

                            var drawAxesLinesAndTicksForSameOrdDimGather = function() {

                                svg.selectAll(".axis").remove();

                                drawXAxisLinesAndTicksForSameOrdDimGather();
                                drawYAxisLinesAndTicksForSameOrdDimGather();
                            }

                            var drawXAxisLinesAndTicksForScatter = function() {

                                xAxis = d3.svg.axis()
                                    .scale(xScale)
                                    .ticks(tickGenerator(scope.xdim))
                                    .tickFormat(labelGenerator(scope.xdim))
                                    .orient("bottom");


                                xAxisNodes = svgGroup.append("g")
                                    .attr("class", "x axis")
                                    .attr("transform", "translate(0," + (height) + ")")
                                    .call(xAxis);

                                xAxisNodes.selectAll('text')
                                    .style("font-size", 12);


                                svg.selectAll(".x .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");
                            };

                            var drawYAxisLinesAndTicksForScatter = function() {

                                yAxis = d3.svg.axis()
                                    .scale(yScale)
                                    .ticks(tickGenerator(scope.ydim))
                                    .tickFormat(labelGenerator(scope.ydim))
                                    .orient("left");

                                yAxisNodes = svgGroup.append("g")
                                    .attr("class", "y axis")
                                    .call(yAxis);

                                yAxisNodes.selectAll('text')
                                    .style("font-size", 12);

                                svg.selectAll(".y .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");

                            };




                            var drawXAxisLinesAndTicksForOrdinalGather = function() {

                                var ticks = tickValueGeneratorForOrdinalGather(scope.xdim);

                                xAxis = d3.svg.axis()
                                    .scale(xScale)
                                    .tickValues(ticks)
                                    .tickFormat(labelGeneratorForOrdinalGather(scope.xdim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("bottom");

                                xAxisNodes = svgGroup.append("g")
                                    .attr("class", "x axis")
                                    .attr("transform", "translate(0," + (height) + ")")
                                    .call(xAxis);

                                xAxisNodes.selectAll('text')
                                    .style("font-size", 12);

                                svg.selectAll(".x .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");

                            };

                            var drawYAxisLinesAndTicksForOrdinalGather = function() {

                                var ticks = tickValueGeneratorForOrdinalGather(scope.ydim);

                                yAxis = d3.svg.axis()
                                    .scale(yScale)
                                    .tickValues(ticks)
                                    .tickFormat(labelGeneratorForOrdinalGather(scope.ydim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("left");

                                yAxisNodes = svgGroup.append("g")
                                    .attr("class", "y axis")
                                    .call(yAxis);

                                yAxisNodes.selectAll('text')
                                    .style("font-size", 12);

                                svg.selectAll(".y .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");

                            };

                            var drawXAxisLinesAndTicksForSameOrdDimGather = function() {

                                var ticks = tickValueGeneratorForOrdinalGather(scope.xdim);

                                var calculatedPositions = getCalculatedPositions(scope.xdim);

                                var domain = [calculatedPositions[0], calculatedPositions[calculatedPositions.length - 1]];


                                var xScaleForSameOrdDimGather = d3.scale.linear().domain(domain).range([0, width]);

                                xAxis = d3.svg.axis()
                                    .scale(xScaleForSameOrdDimGather)
                                    .tickValues(ticks)
                                    .tickFormat(labelGeneratorForOrdinalGather(scope.xdim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("bottom");

                                xAxisNodes = svgGroup.append("g")
                                    .attr("class", "x axis")
                                    .attr("transform", "translate(0," + (height) + ")")
                                    .call(xAxis);

                                xAxisNodes.selectAll('text')
                                    .style("font-size", 12);

                                svg.selectAll(".x .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");

                            };

                            var drawYAxisLinesAndTicksForSameOrdDimGather = function() {


                                var ticks = tickValueGeneratorForOrdinalGather(scope.ydim);

                                var calculatedPositions = getCalculatedPositions(scope.xdim);

                                var domain = [calculatedPositions[0], calculatedPositions[calculatedPositions.length - 1]];


                                var yScaleForSameOrdDimGather = d3.scale.linear().domain(domain).range([height, 0])

                                yAxis = d3.svg.axis()
                                    .scale(yScaleForSameOrdDimGather)
                                    .tickValues(ticks)
                                    .tickFormat(labelGeneratorForOrdinalGather(scope.ydim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("left");

                                yAxisNodes = svgGroup.append("g")
                                    .attr("class", "y axis")
                                    .call(yAxis);

                                yAxisNodes.selectAll('text')
                                    .style("font-size", 12);

                                svg.selectAll(".y .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "black");


                            };



                            //returns path string d for <path d="This string">
                            //a curly brace between x1,y1 and x2,y2, w pixels wide 
                            //and q factor, .5 is normal, higher q = more expressive bracket 
                            var makeCurlyBrace = function(x1, y1, x2, y2, w, q) {
                                //Calculate unit vector
                                var dx = x1 - x2;
                                var dy = y1 - y2;
                                var len = Math.sqrt(dx * dx + dy * dy);

                                if (len === 0) {
                                    dx = 0;
                                    dy = 0;
                                } else {
                                    dx = dx / len;
                                    dy = dy / len;
                                }
                                //Calculate Control Points of path,
                                var qx1 = x1 + q * w * dy;
                                var qy1 = y1 - q * w * dx;
                                var qx2 = (x1 - .25 * len * dx) + (1 - q) * w * dy;
                                var qy2 = (y1 - .25 * len * dy) - (1 - q) * w * dx;
                                var tx1 = (x1 - .5 * len * dx) + w * dy;
                                var ty1 = (y1 - .5 * len * dy) - w * dx;
                                var qx3 = x2 + q * w * dy;
                                var qy3 = y2 - q * w * dx;
                                var qx4 = (x1 - .75 * len * dx) + (1 - q) * w * dy;
                                var qy4 = (y1 - .75 * len * dy) - (1 - q) * w * dx;

                                return ("M " + x1 + " " + y1 +
                                    " Q " + qx1 + " " + qy1 + " " + qx2 + " " + qy2 +
                                    " T " + tx1 + " " + ty1 +
                                    " M " + x2 + " " + y2 +
                                    " Q " + qx3 + " " + qy3 + " " + qx4 + " " + qy4 +
                                    " T " + tx1 + " " + ty1);
                            };

                            var drawAxesLinesAndTicksForGather = function() {

                                svg.selectAll(".axis").remove();

                                drawXAxisLinesAndTicksForGather();
                                drawYAxisLinesAndTicksForGather();

                            };

                            var drawXAxisLinesAndTicksForGather = function() {

                                if (getDimType(scope.xdim) !== 'ordinal' || findTypeOfXYDim() === 'OrdOrd') {

                                    drawXAxisLinesAndTicksForNominalGather();
                                } else {

                                    drawXAxisLinesAndTicksForOrdinalGather();
                                }

                            };

                            var drawYAxisLinesAndTicksForGather = function() {



                                if (getDimType(scope.ydim) !== 'ordinal' || findTypeOfXYDim() === 'OrdOrd') {

                                    drawYAxisLinesAndTicksForNominalGather();
                                } else {

                                    drawYAxisLinesAndTicksForOrdinalGather();
                                }


                            };

                            var drawXAxisLinesAndTicksForNominalGather = function() {

                                xAxis = d3.svg.axis()
                                    .scale(xScale)
                                    .tickValues(tickValueGeneratorForGather(scope.xdim))
                                    .tickFormat(labelGeneratorForGather(scope.xdim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("bottom");

                                svg.selectAll(".axis").remove();

                                xAxisNodes = svgGroup.append("g")
                                    .attr("class", "x axis")
                                    .attr("transform", "translate(0," + (height) + ")")
                                    .call(xAxis);

                                xAxisNodes.selectAll('text')
                                    .style("font-size", 10);

                                d3.selectAll(".x .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "white");

                                var xAxisBracketGroup = xAxisNodes.selectAll(".tick")
                                    .append("g")
                                    .attr("x", xBracketGroup)
                                    .attr("y", 0)
                                    .attr("class", "x controlButtonBracketGroup")
                                    .attr("width", widthBracketGroup)
                                    .attr("height", 30)
                                    .attr("rx", 5)
                                    .attr("ry", 5);

                                if (scope.config.isInteractiveAxis) {



                                    xAxisBracketGroup
                                        .on("mouseover", function(d) {
                                            d3.select(this).selectAll("rect")
                                                .style("opacity", 0.7);
                                            d3.select(this).selectAll("text")
                                                .style("opacity", 0.7);
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).selectAll("rect")
                                                .transition()
                                                .duration(1500)
                                                .style("opacity", 0);

                                            d3.select(this).selectAll("text")
                                                .transition()
                                                .duration(1500)
                                                .style("opacity", 0);
                                        });



                                    xAxisBracketGroup.append("text")
                                        .style("opacity", 0)
                                        .style("fill", "black")
                                        .attr("x", 0)
                                        .attr("y", 60 - 30)
                                        .attr("class", "x controlButtonBracket")
                                        .attr("width", widthBracketGroup)
                                        .attr("height", 10)
                                        .attr("dy", 10)
                                        .style("text-anchor", "middle")
                                        .text("Minimize");

                                    xAxisBracketGroup.append("text")
                                        .style("opacity", 0)
                                        .style("fill", "black")
                                        .attr("x", 0)
                                        .attr("y", 60 - 14)
                                        .attr("class", "x controlButtonBracket")
                                        .attr("width", widthBracketGroup)
                                        .attr("height", 10)
                                        .attr("dy", 10)
                                        .style("text-anchor", "middle")
                                        .text("Maximize");


                                    //     });

                                    xAxisBracketGroup.append("rect")
                                        .style("opacity", 0)
                                        .style("fill", "gray")
                                        .attr("x", xBracketGroup)
                                        .attr("y", 60 - 32)
                                        .attr("class", "x controlButtonBracket")
                                        .attr("width", widthBracketGroup)
                                        .attr("height", 14)
                                        .attr("rx", 5)
                                        .attr("ry", 5)
                                        .on("mouseover", function(d) {
                                            d3.select(this).style("fill", 'lightsteelblue');
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).style("fill", 'lightgray')

                                        })
                                        .on("click", function(d, i) {

                                            toggleMinimizeCluster(scope.xdim, i);
                                        });

                                    xAxisBracketGroup.append("rect")
                                        .style("opacity", 0)
                                        .style("fill", "gray")
                                        .attr("x", xBracketGroup)
                                        .attr("y", 60 - 16)
                                        .attr("class", "x controlButtonBracket")
                                        .attr("width", widthBracketGroup)
                                        .attr("height", 14)
                                        .attr("rx", 5)
                                        .attr("ry", 5)
                                        .on("mouseover", function(d) {
                                            d3.select(this).style("fill", 'green');
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).style("fill", 'lightgray')

                                        })
                                        .on("click", function(d, i) {
                                            //console.log(d);
                                            // toggleMinimizeCluster(scope.xdim, i);
                                            toggleMaximizeCluster(scope.xdim, i)
                                        });



                                }




                                xAxisBracketGroup.append("path")
                                    .attr("class", "x bracket")
                                    .transition()
                                    .duration(500)
                                    .attr("d", pathXBracket);




                            };


                            var drawYAxisLinesAndTicksForNominalGather = function() {




                                yAxis = d3.svg.axis()
                                    .scale(yScale)
                                    .tickValues(tickValueGeneratorForGather(scope.ydim))
                                    .tickFormat(labelGeneratorForGather(scope.ydim))
                                    .tickSize(12, 0) //Provides 0 size ticks at center position for gather
                                    .orient("left");


                                yAxisNodes = svgGroup.append("g")
                                    .attr("class", "y axis")
                                    .call(yAxis);


                                yAxisNodes.selectAll('text')
                                    .style("font-size", 10);

                                d3.selectAll(".y .tick line")
                                    .style("stroke-width", 1)
                                    .style("stroke", "white");


                                var yAxisBracketGroup = yAxisNodes.selectAll(".tick")
                                    .append("g")
                                    .attr("x", 0)
                                    .attr("y", yBracketGroup)
                                    .attr("class", "y controlButtonBracketGroup")
                                    .attr("width", margin)
                                    .attr("height", heightBracketGroup)
                                    .attr("rx", 5)
                                    .attr("ry", 5);



                                if (scope.config.isInteractiveAxis) {

                                    yAxisBracketGroup
                                        .on("mouseover", function(d) {
                                            d3.select(this).selectAll("rect")
                                                .style("opacity", 0.9);
                                            d3.select(this).selectAll("text")
                                                .style("opacity", 0.9);
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).selectAll("rect")
                                                .transition()
                                                .duration(2000)
                                                .style("opacity", 0);

                                            d3.select(this).selectAll("text")
                                                .transition()
                                                .duration(2000)
                                                .style("opacity", 0);
                                        });



                                    yAxisBracketGroup.append("text")
                                        .style("opacity", 0)
                                        .style("fill", "black")
                                        .attr("x", 20)
                                        .attr("y", 0)
                                        .attr("class", "y controlButtonBracket")
                                        .attr("width", 20)
                                        .attr("height", heightBracketGroup)
                                        .attr("dy", 10)
                                        .style("text-anchor", "left")
                                        .text("Minimize");

                                    yAxisBracketGroup.append("text")
                                        .style("opacity", 0)
                                        .style("fill", "black")
                                        .attr("x", 110)
                                        .attr("y", 0)
                                        .attr("class", "y controlButtonBracket")
                                        .attr("width", 10)
                                        .attr("height", heightBracketGroup)
                                        .attr("dy", 10)
                                        .style("text-anchor", "left")
                                        .text("Maximize");


                                    //     });

                                    yAxisBracketGroup.append("rect")
                                        .style("opacity", 0)
                                        .style("fill", "gray")
                                        .attr("x", 10)
                                        .attr("y", -2)
                                        .attr("class", "y controlButtonBracket")
                                        .attr("width", margin)
                                        .attr("height", 14)
                                        .attr("rx", 5)
                                        .attr("ry", 5)
                                        .on("mouseover", function(d) {
                                            d3.select(this).style("fill", 'lightsteelblue');
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).style("fill", 'lightgray')

                                        })
                                        .on("click", function(d, i) {

                                            toggleMinimizeCluster(scope.ydim, i);
                                        });

                                    yAxisBracketGroup.append("rect")
                                        .style("opacity", 0)
                                        .style("fill", "gray")
                                        .attr("x", 100)
                                        .attr("y", -2)
                                        .attr("class", "y controlButtonBracket")
                                        .attr("width", margin)
                                        .attr("height", 14)
                                        .attr("rx", 5)
                                        .attr("ry", 5)
                                        .on("mouseover", function(d) {
                                            d3.select(this).style("fill", 'green');
                                        })
                                        .on("mouseout", function(d) {


                                            d3.select(this).style("fill", 'lightgray')

                                        })
                                        .on("click", function(d, i) {
                                            console.log(d);
                                            // toggleMinimizeCluster(scope.xdim, i);
                                            toggleMaximizeCluster(scope.ydim, i)
                                        });

                                }

                                yAxisNodes.selectAll(".tick")
                                    .append("path")
                                    .attr("class", "y bracket")
                                    .transition()
                                    .duration(500)
                                    .attr("d", pathYBracket);



                            };

                            var toggleMinimizeCluster = function(dim, i) {


                                var key = d3.map(dimSetting[dim].keyValue).values()[i].keyValue;

                                var keyObject = dimSetting[dim].keyValue[key];

                                keyObject.isMinimized = !keyObject.isMinimized;

                                drawPlot();

                            };

                            var toggleMaximizeCluster = function(dim, i) {


                                var key = d3.map(dimSetting[dim].keyValue).values()[i].keyValue;

                                var keyObject = dimSetting[dim].keyValue[key];

                                keyObject.isMaximized = !keyObject.isMaximized;

                                var keyValue = d3.map(dimSetting[dim].keyValue).values();


                                if (keyObject.isMaximized === true) {


                                    keyValue.forEach(function(d) {

                                        d.isMinimized = true;


                                    });

                                    keyObject.isMinimized = false;


                                } else {
                                    keyValue.forEach(function(d) {

                                        d.isMinimized = false;


                                    });

                                }

                                drawPlot();

                            };

                            var pathXBracket = function(d, i) {

                                var dim = scope.xdim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfCluster(dim, key, xScale);

                                if (length === 0) {
                                    return ("M 0 0 " +
                                        " L 0 " + 10);
                                } else {

                                    return makeCurlyBrace(-length / 2, 2, length / 2, 2, 10, 0.6);
                                }
                            };

                            var pathYBracket = function(d, i) {

                                var dim = scope.ydim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfCluster(dim, key, yScale);

                                if (length === 0) {
                                    return ("M 0 0 " +
                                        " L -10 " + 0);
                                } else {

                                    return makeCurlyBrace(-2, length / 2, -2, -length / 2, 10, 0.6);
                                }



                            };


                            var xBracket = function(d, i) {

                                var dim = scope.xdim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfCluster(dim, key, xScale);

                                return length / 2 * (-1);

                            };

                            var xBracketGroup = function(d, i) {

                                var dim = scope.xdim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfClusterIncludingMargin(dim, key, xScale);

                                return length / 2 * (-1);

                            };

                            var widthBracket = function(d, i) {

                                var dim = scope.xdim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfCluster(dim, key, xScale);

                                return length;

                            };

                            var widthBracketGroup = function(d, i) {

                                var dim = scope.xdim;

                                var key = getKeyFromIndex(dim, i);

                                var length = lengthOfClusterIncludingMargin(dim, key, xScale);

                                return length;

                            };

                            var yBracket = function(d, i) {

                                var dim = scope.ydim;

                                var key = getKeyFromIndex(dim, i);

                                var length = -lengthOfCluster(dim, key, yScale);

                                return length / 2 * (-1);

                            };

                            var yBracketGroup = function(d, i) {

                                var dim = scope.ydim;

                                var key = getKeyFromIndex(dim, i);

                                var length = -lengthOfClusterIncludingMargin(dim, key, yScale);

                                return length / 2 * (-1);

                            };

                            var heightBracket = function(d, i) {

                                var dim = scope.ydim;

                                var key = getKeyFromIndex(dim, i);

                                var length = -lengthOfCluster(dim, key, yScale);

                                return length;

                            };

                            var heightBracketGroup = function(d, i) {

                                var dim = scope.ydim;

                                var key = getKeyFromIndex(dim, i);

                                var length = -lengthOfClusterIncludingMargin(dim, key, yScale);

                                return length;

                            };


                            var lengthOfCluster = function(dim, key, scale) {

                                var keyObject = dimSetting[dim].keyValue[key];

                                if (keyObject.isMinimized) {

                                    return 0;

                                } else {

                                    return scale(1 - 2 * marginClusterRatio) - scale(0);
                                }



                            };

                            var lengthOfClusterIncludingMargin = function(dim, key, scale) {

                                var keyObject = dimSetting[dim].keyValue[key];

                                if (keyObject.isMinimized) {

                                    return scale(2 * marginClusterRatio) - scale(0);

                                } else {

                                    return scale(1) - scale(0);
                                }



                            };



                            var getKeyFromIndex = function(dim, i) {

                                if (!dimSetting[dim].keyValue) {

                                    debugger;
                                    console.log(dim);
                                }
                                if (!d3.map(dimSetting[dim].keyValue).values()[i]) {

                                    debugger;
                                    console.log(dim);
                                }

                                return d3.map(dimSetting[dim].keyValue).values()[i].keyValue;

                            };


                            var drawAxesLabel = function() {

                                xAxisNodes
                                    .append("text")
                                    .attr("class", "axislabel")
                                    .attr("x", width / 2)
                                    .attr("y", 56)
                                    .style("text-anchor", "middle")
                                    .text(scope.xdim);

                                //Setup Y axis

                                yAxisNodes
                                    .append("text")
                                    .attr("class", "axislabel")
                                    .style("text-anchor", "middle")
                                    .attr('transform', function(d, i) { // NEW
                                        var vert = height / 2; // NEW
                                        // var horz = -margin / 2; // NEW
                                        var horz = -60;
                                        return 'translate(' + horz + ',' + vert + ')rotate(-90)'; // NEW
                                    })
                                    .text(scope.ydim);


                                // yAxisNodes
                                //     .append("text")
                                //     .attr("class", "axislabel")
                                //     .text(findDisplayName(scope.ydim))
                                //     .attr('transform', function(d, i) { // NEW
                                //         var vert = height / 2; // NEW
                                //         var horz = -margin / 2; // NEW
                                //         return 'translate(' + horz + ',' + vert + ')rotate(-90)'; // NEW
                                //     });



                            };

                            var drawLegends = function() {

                                resetLegends();

                                if (!scope.config.colorDim) {

                                    return;
                                }

                                var currentDimSetting = dimSetting[scope.config.colorDim];

                                if (currentDimSetting.dimType === 'ordinal') {

                                    drawHeatMapLegends();
                                } else {

                                    drawNominalLegends();
                                }
                            };

                            var resetLegends = function() {

                                var legendGroup = svg.selectAll(".legend").remove();

                            };

                            var drawHeatMapLegends = function() {

                                var colorDomain = d3.extent(scope.data, function(d) {
                                    return +d[scope.config.colorDim];
                                });

                                var widthHeatMap = 200;
                                var heightHeatMap = 18;


                                var xScaleForHeatMap = d3.scale.linear()
                                    .domain(colorDomain)
                                    .rangeRound([width - 100, width + 100]);

                                var values = d3.range(colorDomain[0], colorDomain[1], (colorDomain[1] - colorDomain[0]) / widthHeatMap);

                                var g = svg.append("g")
                                    .attr("class", "legend");



                                var heatmap = g.selectAll("rect")
                                    .data(values)
                                    .enter().append("rect")
                                    .attr("x", xScaleForHeatMap)
                                    .attr("y", 20)
                                    .attr("width", 1)
                                    .attr("height", heightHeatMap)
                                    .style("fill", colorScaleForHeatMap);

                                g.append("text")
                                    .attr("x", width + 12)
                                    .attr("y", 10)
                                    .attr("dy", ".35em")
                                    .style("text-anchor", "middle")
                                    .text(scope.config.colorDim);

                                g.append("text")
                                    .attr("x", xScaleForHeatMap(values[0]))
                                    .attr("y", 50)
                                    .attr("dy", ".35em")
                                    .style("text-anchor", "middle")
                                    .text(d3.round(colorDomain[0], 1));

                                g.append("text")
                                    .attr("x", xScaleForHeatMap(values[values.length - 1]))
                                    .attr("y", 50)
                                    .attr("dy", ".35em")
                                    .style("text-anchor", "middle")
                                    .text(d3.round(colorDomain[1], 1));

                            };

                            var drawNominalLegends = function() {


                                var legendGroup = svg.selectAll(".legend")
                                    .data(getKeys(scope.config.colorDim), function(d) {
                                        return d;
                                    });

                                legendGroup.exit().remove();


                                var legend = legendGroup.enter().append("g")
                                    .attr("class", "legend")
                                    .attr("transform", function(d, i) {
                                        return "translate(0," + (i * 20 + 5) + ")";
                                    });

                                legend.append("rect")
                                    .attr("x", width - 18)
                                    .attr("width", 18)
                                    .attr("height", 18)
                                    .style("fill", function(d) {
                                        return color(d);
                                    });

                                legend.append("text")
                                    .attr("x", width + 5)
                                    .attr("y", 9)
                                    .attr("dy", ".35em")
                                    .style("text-anchor", "left")
                                    .text(function(d) {
                                        return d;
                                    });



                                var g = svg.append("g")
                                    .attr("class", "legend");



                                g.append("text")
                                    .attr("x", width - 24)
                                    .attr("y", 10)
                                    .attr("dy", ".35em")
                                    .style("text-anchor", "end")
                                    .text(scope.config.colorDim);




                            }; //End renderer

                        }

                    }; //End return 

                }
            ] // End function (d3Service)

        );

    angular.module('myApp.directives')
        .directive('focusMe', function($timeout, $parse) {
            return {
                //scope: true,   // optionally create a child scope
                link: function(scope, element, attrs) {
                    var model = $parse(attrs.focusMe);
                    scope.$watch(model, function(value) {
                        console.log('value=', value);
                        if (value === true) {
                            $timeout(function() {
                                element[0].focus();
                            });
                        }
                    });
                    // to address @blesh's comment, set attribute value to 'false'
                    // on blur event:
                    element.bind('blur', function() {
                        console.log('blur');
                        scope.$apply(model.assign(scope, false));
                    });
                }
            };
        });

}());
