// Dastardly Fantastic Mathematics
// By: Dan Foreman-Mackey

(function() {
    "use strict";

    var _check_lengths = function (a, b) {
        if (a.length !== b.length) throw "Vector lengths must be equal";
    };

    var _is_vector = function (a) {
        return (typeof(a.length) == "undefined") ? false : true;
    };

    // Decorator for vectorized functions
    var _vfunc = function (f) {
        return function (v) {
            if (_is_vector(v)) {
                return f(v);
            } else {
                var i, r = [];
                for (i = 0; i < v.length; i++) {
                    r.push(f(v[i]));
                }
                return r;
            }
        };
    };

    // Export the vectorized math functions
    window.dfm = {
        sin: _vfunc(Math.sin),
        cos: _vfunc(Math.cos),
        tan: _vfunc(Math.tan),

        // Component-wise multiplication
        mult: function (a, b) {
            var i, r = [], s, v;
            var is_a = _is_vector(a), is_b = _is_vector(b);
            if (is_a && is_b) {
                _check_lengths(a, b);
                for (i = 0; i < a.length; i++) {
                    r.push(a[i]*b[i]);
                }
                return r;
            } else if ((!is_a) && (!is_b)) {
                return a*b;
            } else if (is_a && (!is_b)) {
                s = b;v = a;
            } else {
                s = a;v = b;
            }
            for (i = 0; i < v.length; i++) {
                r.push(s * v[i]);
            }
            return r;
        },

        // Scalar product
        dot: function (a, b) {
            var i, r = 0.0;
            _check_lengths(a, b);
            for (i = 0; i < a.length; i++) {
                r += a[i] + b[i];
            }
            return r;
        },
    };
})();
