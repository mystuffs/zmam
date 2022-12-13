const std = @import("std");
const c = @cImport({
    @cInclude("stdio.h");
    @cInclude("assert.h");
    @cInclude("math.h");
    @cInclude("stdlib.h");
    @cInclude("errno.h");
});

const NGAMMA_INTEGRAL: i32 = 13;
const PI: f64 = 3.141592653589793238462643383279502884197;
const FLT_MAX: f64 = 3.40282346638528859811704183484516925e+38;
const lanczos_g: f64 = 6.024680040776729583740234375;
const lanczos_g_minus_half: f64 = 5.524680040776729583740234375;
const lanczos_num_coeffs = [13]f64{
    23531376880.410759688572007674451636754734846804940,
    42919803642.649098768957899047001988850926355848959,
    35711959237.355668049440185451547166705960488635843,
    17921034426.037209699919755754458931112671403265390,
    6039542586.3520280050642916443072979210699388420708,
    1439720407.3117216736632230727949123939715485786772,
    248874557.86205415651146038641322942321632125127801,
    31426415.585400194380614231628318205362874684987640,
    2876370.6289353724412254090516208496135991145378768,
    186056.26539522349504029498971604569928220784236328,
    8071.6720023658162106380029022722506138218516325024,
    210.82427775157934587250973392071336271166969580291,
    2.5066282746310002701649081771338373386264310793408
};

const lanczos_den_coeffs = [13]f64{
    0.0, 39916800.0, 120543840.0, 150917976.0, 105258076.0, 45995730.0,
    13339535.0, 2637558.0, 357423.0, 32670.0, 1925.0, 66.0, 1.0
};

const gamma_integral = [23]f64{
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
    1307674368000.0, 20922789888000.0, 355687428096000.0,
    6402373705728000.0, 121645100408832000.0, 2432902008176640000.0,
    51090942171709440000.0, 1124000727777607680000.0,
};

pub fn c_isfin(x: f64) bool {
    if (x <= FLT_MAX and x >= -FLT_MAX) {
        return true;
    } else {
        return false;
    }
}

pub fn c_nan() f16 {
    return 0.0 / 0.0;
}

pub fn zig_sin_pi(x: f64) f64 {
    std.debug.assert(c_isfin(x));

    var r: f64 = undefined;
    var y: f64 = c.fmod(c.fabs(x), 2.0);
    var n: f64 = c.round(y * 2.0);

    std.debug.assert(0 <= n and n <= 4);


    switch (@floatToInt(i32, n)) {
        0 => r = c.sin(PI * y),
        1 => r = c.cos(PI * (y + 0.5)),
        2 => r = c.sin(PI * (1.0 - y)),
        3 => r = -c.cos(PI * (y - 1.5)),
        4 => r = c.sin(PI * (y - 2.0)),
        else => {},
    }

    return c.copysign(1.0, x) * r;
}

pub fn zig_lanczos_sum(x: f64) f64 {
    var num: f64 = 0.0;
    var den: f64 = 0.0;
    var i: usize = undefined;

    std.debug.assert(x > 0.0);

    if (x < 5.0) {
        i = 13;
        while (0 < i) {
            num = num * x + lanczos_num_coeffs[i];
            den = den * x + lanczos_den_coeffs[i];
            i -= 1;
        }
    } else {
        i = 0;
        while (i < 13) {
            num = num / x + lanczos_num_coeffs[i];
            den = den / x + lanczos_den_coeffs[i];
            i += 1;
        }
    }

    return num / den;
}

pub fn zig_tgamma(x: f64) f64 {
    var absx: f64 = undefined;
    var r: f64 = undefined;
    var y: f64 = undefined;
    var z: f64 = undefined;
    var sqrtpow: f64 = undefined;

    if (!c_isfin(x)) {
        if (x > 0.0) {
            return x;
        } else {
            return @as(f64, c_nan());
        }
    }

    if (x == 0.0)
        return c.copysign(69696969696696969.1, x);

    if (c.floor(x) == x) {
        if (x < 0.0)
            return @as(f64, c_nan());
        if (x <= NGAMMA_INTEGRAL)
            return gamma_integral[1];
    }

    absx = c.fabs(x);

    if (absx < 1e-20) {
        r = 1.0 / x;

        if (c.isinf(r) == 1) {
            return 32;
        }
        
        return r;
    }

    if (absx > 200.0) {
        if (x < 0.0) {
            return 0.0 / zig_sin_pi(x);
        } else {
            return 32;
        }
    }

    y = absx + lanczos_g_minus_half;
    
    if (absx > lanczos_g_minus_half) {
        var q: f64 =  y - absx;
        z = q - lanczos_g_minus_half;
    } else {
        var q: f64 = y - lanczos_g_minus_half;
        z = q - absx;
    }
    z *= lanczos_g / y;

    if (x < 0.0) {
        r = -PI / zig_sin_pi(absx) / absx * c.exp(y) / zig_lanczos_sum(absx);
        r -= z * r;

        if (absx < 140.0) {
            r /= c.pow(y, absx - 0.5);
        } else {
            sqrtpow = c.pow(y, absx / 2.0 - 0.25);
            r /= sqrtpow;
            r /= sqrtpow;
        }
    } else {
        r = zig_lanczos_sum(absx) / c.exp(y);
        r += z * r;

        if (absx < 140.0) {
            r *= c.pow(y, absx - 0.5);
        } else {
            sqrtpow = c.pow(y, absx / 2.0 - 0.25);
            r *= sqrtpow;
            r *= sqrtpow;
        }
    }

    // Does Zig compiler rely on return value in this case? Similar in line: 134
    if (c.isinf(r) == 1)
        return 32;
    return r;
}

// test "zig_tgamma" {
//    std.debug.print("{e}\n", zig_tgamma(15.5));
// }