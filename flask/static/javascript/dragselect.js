var selectionRect = {
	element			: null,
	previousElement : null,
	currentY		: 0,
	currentX		: 0,
	originX			: 0,
	originY			: 0,
	setElement: function(ele) {
		this.previousElement = this.element;
		this.element = ele;
	},
	getNewAttributes: function() {
		var x = this.currentX<this.originX?this.currentX:this.originX;
		var y = this.currentY<this.originY?this.currentY:this.originY;
		var width = Math.abs(this.currentX - this.originX);
		var height = Math.abs(this.currentY - this.originY);
		return {
			x       : x,
			y       : y,
			width  	: width,
			height  : height
		};
	},
	getCurrentAttributes: function() {
		// use plus sign to convert string into number
		var x = +this.element.attr("x");
		var y = +this.element.attr("y");
		var width = +this.element.attr("width");
		var height = +this.element.attr("height");
		return {
			x1  : x,
			y1	: y,
			x2  : x + width,
			y2  : y + height
		};
	},
	init: function(newX, newY) {
		var rectElement = svg.append("rect")
			.attr({
				rx      : 4,
				ry      : 4,
				x       : 0,
				y       : 0,
				width   : 0,
				height  : 0
			})
			.classed("selection", true);
		this.setElement(rectElement);
		this.originX = newX;
		this.originY = newY;
		this.update(newX, newY);
	},
	update: function(newX, newY) {
		this.currentX = newX;
		this.currentY = newY;
		this.element.attr(this.getNewAttributes());
	},
	focus: function() {
		this.element
			.style("stroke", "#DE695B")
			.style('fill-opacity', 0)
			.style("stroke-width", "1");
	},
	remove: function() {
		this.element.remove();
		this.element = null;
	},
	removePrevious: function() {
		if(this.previousElement) {
			this.previousElement.remove();
		}
	}
};