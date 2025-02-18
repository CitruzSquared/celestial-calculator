var pi = Math.PI;
var G = 6.674e-20;

var earth_sma = 149598023.0;
var earth_radius = 6378.137;
var earth_mass = 5.972168e24;
var earth_eccentricity = 0.0167086;
var earth_AOP = 1.79676742118;
var day = 24.0;

var sun_radius = 695700.0;
var sun_mass = 1.9885e30;

var moon_sma = 384399.0;
var moon_radius = 1737.4;
var moon_mass = 7.346e22;
var moon_eccentricity = 0.0549;
var moon_inclination = 0.08979719001;

var body = document.body;



var earth_period = period(sun_mass, earth_mass, earth_sma);
var moon_period = period(earth_mass, moon_mass, moon_sma);
var moon_synodic = earth_period * moon_period / (earth_period - moon_period);
var earth_periapsis = earth_sma * (1 - earth_eccentricity);
var earth_apoapsis = earth_sma * (1 + earth_eccentricity);
var moon_periapsis = moon_sma * (1 - moon_eccentricity);
var moon_apoapsis = moon_sma * (1 + moon_eccentricity);

body.innerHTML += "<h2> Basic Orbital Properties </h2>";
body.innerHTML += "<p> Earth periapsis: " + earth_periapsis.toString() + " km </p>";
body.innerHTML += "<p> Earth apoapsis: " + earth_apoapsis.toString() + " km </p>";
body.innerHTML += "<p> Moon periapsis: " + moon_periapsis.toString() + " km </p>";
body.innerHTML += "<p> Moon apoapsis: " + moon_apoapsis.toString() + " km </p>";
body.innerHTML += "<p> Earth period: " + earth_period.toString() + " dy </p>";
body.innerHTML += "<p> Moon sidereal period: " + moon_period.toString() + " dy </p>";
body.innerHTML += "<p> Moon synodic period: " + moon_synodic.toString() + " dy </p>";

body.innerHTML += "<h2> Calendar Properties </h2>";
body.innerHTML += "<p> Solar calendar leaps: " + dec_to_frac(earth_period).join(", ") + "</p>";
body.innerHTML += "<p> Lunar calendar leaps: " + dec_to_frac(moon_synodic).join(", ") + "</p>";
body.innerHTML += "<p> Lunisolar calendar leaps: " + dec_to_frac(earth_period/moon_synodic).join(", ") + "</p>";

function period(mass1, mass2, sma) {
	return Math.sqrt(4 * pi * pi * sma * sma * sma / (G * (mass1 + mass2))) / 3600 / day;
}

function dec_to_frac(number) {
	let decimal = number - Math.floor(number);
	let denoms = [];
	let fracs = [];
	for (let i = 0; i < 10; i++) {
		let reciprocal = 1 / decimal;
		denoms.push(Math.floor(reciprocal));
		decimal = reciprocal - Math.floor(reciprocal);
	}
	for (let i = 0; i < 10; i++) {
		let numerator = 1;
		let denominator = denoms[i];
		for (let j = i - 1; j >= 0; j--) {
			numerator_copy = numerator;
			numerator = denominator;
			denominator = denoms[j] * denominator + numerator_copy;
		}
		
		fracs.push(numerator.toString() + "/" + denominator.toString());
	}
	return(fracs);
}

console.log(dec_to_frac(pi));
