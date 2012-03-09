(function () {
    "use strict";

    // Some math helpers
    window.radians = function (d) {return d * Math.PI/180.0;};
    window.linspace = function (a, b, N) {
        var i, r = [];
        if (a >= b) throw "Don't be silly";
        for (i = 0; i < N; i++)
            r.push(a + i*(b-a)/(N-1));
        return r;
    };
    window.arange = function (a, b, d) {
        var i, r = [];
        if (!d) d = 1;
        if (typeof(b) == "undefined") b = a, a = 0;
        for (i = a; i < b; i += d)
            r.push(i);
        return r;
    };

    // A single planet object
    var Planet = function (Rp, R, incl) {
        this.Rp = Rp;
        this.R = R;
        this.incl = radians(incl);
        this._si = Math.sin(this.incl);
        this._si2 = this._si * this._si;
        this._ci = Math.cos(this.incl);
        return this;
    };
    Planet.prototype = {
        constructor: Planet,

        position: function (ph) {
            var s = Math.sin(ph), c = Math.cos(ph), proj, depth;
            proj = this.R * Math.sqrt(s*s + c*c*this._si2);
            depth = this.R*this._ci*c;
            return [proj, depth];
        },
    };

    // A full planetary system
    var PlanetarySystem = function (Rs) {
        this.Rs = Rs;
        this.planets = [];
        return this;
    };
    PlanetarySystem.prototype = {
        constructor: PlanetarySystem,

        add_planet: function (p) {
            this.planets.push(p);
        },

        occultation: function (z, p) {
            if (Math.abs(1-p) < z && z < 1+p) {
                var p2 = p*p, z2 = z*z, k1, k2, k3;
                var tmp = z2-p2+1;

                k1 = Math.acos(0.5*tmp/z);
                k2 = Math.acos(0.5*(z2+p2-1)/z/p);
                k3 = Math.sqrt(z2-0.25*tmp*tmp);

                return (k1 + p2 * k2 - k3)/Math.PI;
            } else if (z <= 1-p) return p*p;
            else if (z <= p-1) return 1;
            return 0;
        },

        lightcurve: function (t) {
            var me = this;
            return _.map(t, function (t0) {
                return {t: t0, f: 1-_.reduce(me.planets, function (v, p) {
                        var pos, proj;
                        console.log(p);
                        pos = p.position(t0);
                        proj = pos[0]/me.Rs;
                        return v+me.occultation(proj, p.Rp/me.Rs);
                    }, 0)
                };
            });
        },
    };

    window.Planet = Planet;
    window.PlanetarySystem = PlanetarySystem;
})();

$(function () {
    var m = [80, 80, 80, 80],
        w = 960 - m[1] - m[3],
        h = 500 - m[0] - m[2];

    // Scales and axes. Note the inverted domain for the y-scale: bigger is up!
    var x = d3.scale.linear().range([0, w]),
        y = d3.scale.linear().range([h, 0]),
        xAxis = d3.svg.axis().scale(x).ticks(5),
        yAxis = d3.svg.axis().scale(y).ticks(4).orient("left");

    // A line generator, for the dark stroke.
    var line = d3.svg.line()
        .x(function(d) { return x(d.t); })
        .y(function(d) { return y(d.f); });

    var p = new Planet(0.6, 10, 6), ps = new PlanetarySystem(1);
    ps.add_planet(p);
    var lc = ps.lightcurve(arange(-Math.PI, Math.PI, 0.01));

    x.domain([-1, 1]);
    y.domain([_.min(lc, function (v) {return v.f;}).f-0.1, 1.1]).nice();

    var svg = d3.select("#data").append("svg:svg")
        .attr("width", w + m[1] + m[3])
        .attr("height", h + m[0] + m[2])
        .append("svg:g")
        .attr("transform", "translate(" + m[3] + "," + m[0] + ")");

    // Add the clip path.
    svg.append("svg:clipPath")
        .attr("id", "clip")
        .append("svg:rect")
        .attr("width", w)
        .attr("height", h);

    // Add the x-axis.
    svg.append("svg:g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + h + ")")
        .call(xAxis);

    // Add the y-axis.
    svg.append("svg:g")
        .attr("class", "y axis")
        .call(yAxis);

    // Add the line path.
    svg.append("svg:path")
        .attr("class", "line")
        .attr("clip-path", "url(#clip)")
        .attr("d", line(lc));
});
