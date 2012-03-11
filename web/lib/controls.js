(function () {
    "use strict"

    var slider = function (el, callback, range) {
        var this_ = this;
        this.el = el;
        this.callback = callback;
        this.w = 150;
        this.h = 10;
        if (typeof(range) !== "undefined") this.range = range;
        else this.range = [0, 1];
        this._value = 0.5*(range[0]+range[1]);

        this.svg = d3.select(this.el)
            .append("svg:svg")
            .attr("width", this.w)
            .attr("height", this.h)
            .on("mousedown", function () {
                // This is a hacky way of making the mouse events work...
                // not sure why we need it.
                this_.box = this;
                D3SliderController.set_current(this_);
            });
        this.svg.append("svg:clipPath")
            .append("svg:rect")
            .attr("width", this.w)
            .attr("height", this.h);

        this.x = d3.scale.linear().domain(this.range).range([0, this.w]);
        this.y = d3.scale.linear().domain([0, 2]).range([0, this.h]);

        var line = d3.svg.line()
            .x(function(d) { return this_.x(d); })
            .y(function(d) { return this_.y(1); });
        this.g = this.svg.append("g");
        this.g.append("svg:path")
            .attr("d", line(this.range))
            .attr("class", "line");
        this.circ = this.g.append("circle")
            .attr("r", 4)
            .attr("cx", this_.x(this_._value))
            .attr("cy", this_.y(1));
    };
    slider.prototype = {
        constructor: slider,
        value: function (v) {
            if (typeof(v) !== "undefined") {
                this._value = v;
                this.callback(this);
                return this._value;
            } else return this._value;
        },
        draw: function () {
            this.circ.attr("cx", this.x(this._value));
        },
        slide: function () {
            var val = this.x.invert(d3.mouse(this.box)[0]);
            val = Math.max(val, this.range[0]);
            val = Math.min(val, this.range[1]);
            this.value(val);
            this.draw();
        },
    };

    var controller = function () {
        var this_ = this;
        this.current = false;
        d3.select("body").on("mousemove", function () {
            if (this_.current) {
                this_.current.slide();
            }
            return false;
        }).on("mouseup", function () {
            this_.current = false;
        });
    };
    controller.prototype = {
        constructor: controller,
        set_current: function (slider) {
            this.current = slider;
        },
    };

    window.D3Slider = slider;
    window.D3SliderController = new controller;
})();

