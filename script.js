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

// --------------------

var earth_periapsis = earth_sma * (1 - earth_eccentricity);
var earth_apoapsis = earth_sma * (1 + earth_eccentricity);
var moon_periapsis = moon_sma * (1 - moon_eccentricity);
var moon_apoapsis = moon_sma * (1 + moon_eccentricity);

var earth_period = period(sun_mass, earth_mass, earth_sma);
var moon_period = period(earth_mass, moon_mass, moon_sma);
earth_period = 365.256363004;
var moon_period = 27.321661554
var moon_synodic = earth_period * moon_period / (earth_period - moon_period);
var m =  moon_period / earth_period;
var months_in_year = earth_period / moon_synodic;

var seasons = calculate_seasons();

var precessions = calculate_precessions(m);

body.innerHTML += "<h2> Basic Orbital Properties </h2>";
body.innerHTML += "<p> Earth periapsis: " + earth_periapsis.toString() + " km </p>";
body.innerHTML += "<p> Earth apoapsis: " + earth_apoapsis.toString() + " km </p>";
body.innerHTML += "<p> Moon periapsis: " + moon_periapsis.toString() + " km </p>";
body.innerHTML += "<p> Moon apoapsis: " + moon_apoapsis.toString() + " km </p>";
body.innerHTML += "<p> Earth period: " + earth_period.toString() + " dy </p>";
body.innerHTML += "<p> Moon sidereal period: " + moon_period.toString() + " dy </p>";
body.innerHTML += "<p> Moon synodic period: " + moon_synodic.toString() + " dy </p>";
body.innerHTML += "<p> Earth-Moon period ratio: " + m.toString() + "</p>";
body.innerHTML += "<p> Earth-Moon synodic ratio: " + months_in_year.toString() + "</p>";

body.innerHTML += "<h2> Calendar Properties </h2>";
body.innerHTML += "<p> Solar calendar leaps: " + dec_to_frac(earth_period).join(" or ") + " dy/yr </p>";
body.innerHTML += "<p> Lunar calendar leaps: " + dec_to_frac(moon_synodic).join(" or ") + " dy/mn </p>";
body.innerHTML += "<p> Lunisolar calendar leaps: " + dec_to_frac(earth_period/moon_synodic).join(" or ") + " mn/yr</p>";

body.innerHTML += "<h2> Seasonal Properties </h2>";
body.innerHTML += "<p> Northern Lichun / Imbolc: " + seasons[0][0].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Spring Equinox: " + seasons[0][1].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Lixia / Bealtaine: " + seasons[0][2].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Summer Solstice: " + seasons[0][3].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Liqiu / Lunasa: " + seasons[0][4].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Autumn Equinox: " + seasons[0][5].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Lidong / Samhain: " + seasons[0][6].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Winter Solstice: " + seasons[0][7].toString() + " dy after periapsis </p>";
body.innerHTML += "<p> Northern Astronomical Spring length: " + seasons[1][0].toString() + " dy </p>";
body.innerHTML += "<p> Northern Astronomical Summer length: " + seasons[1][1].toString() + " dy </p>";
body.innerHTML += "<p> Northern Astronomical Autumn length: " + seasons[1][2].toString() + " dy </p>";
body.innerHTML += "<p> Northern Astronomical Winter length: " + seasons[1][3].toString() + " dy </p>";
body.innerHTML += "<p> Northern Solar Spring length: " + seasons[2][0].toString() + " dy </p>";
body.innerHTML += "<p> Northern Solar Summer length: " + seasons[2][1].toString() + " dy </p>";
body.innerHTML += "<p> Northern Solar Autumn length: " + seasons[2][2].toString() + " dy </p>";
body.innerHTML += "<p> Northern Solar Winter length: " + seasons[2][3].toString() + " dy </p>";

body.innerHTML += "<h2> Eclipse Properties </h2>";
body.innerHTML += "<p> Lunar nodal precession: " + precessions[0].toString() + " dy/rev </p>";
body.innerHTML += "<p> Lunar apsidal precession: " + precessions[1].toString() + " dy/rev </p>";
body.innerHTML += "<p> Draconic month: " + precessions[2].toString() + " dy </p>";
body.innerHTML += "<p> Anomalistic month: " + precessions[3].toString() + " dy </p>";
body.innerHTML += "<p> Eclipse year: " + precessions[4].toString() + " dy </p>";
body.innerHTML += "<p> Eclipse cycles: " + dec_to_frac(moon_synodic/(precessions[4]/2)).join(" or ") + " half-E.Y./S.M.</p>";

function period(mass1, mass2, sma) {
	return Math.sqrt(4 * pi * pi * sma * sma * sma / (G * (mass1 + mass2))) / 3600 / day;
}

function dec_to_frac(number) {
	let num_approxs = 10;
	let decimal = number - Math.floor(number);
	let denoms = [];
	let fracs = [];
	for (let i = 0; i < num_approxs; i++) {
		let reciprocal = 1 / decimal;
		denoms.push(Math.floor(reciprocal));
		decimal = reciprocal - Math.floor(reciprocal);
	}
	for (let i = 0; i < num_approxs; i++) {
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

function date_from_longitude(longitude) {
	let heliocentric_longitude = longitude - pi;
	if (heliocentric_longitude < 0) {
		heliocentric_longitude = heliocentric_longitude + 2 * pi;
	}
	let true_anomaly = heliocentric_longitude - earth_AOP;
	if (true_anomaly < 0) {
		true_anomaly += 2 * pi;
	}
	let eccentric = Math.acos((earth_eccentricity + Math.cos(true_anomaly)) / (1 + earth_eccentricity * Math.cos(true_anomaly)));
	if (true_anomaly > pi) {
		eccentric = eccentric * -1;
	}
	let mean_anomaly = eccentric - earth_eccentricity * Math.sin(eccentric);
	let date = mean_anomaly * earth_period / (2 * pi);
	if (date < 0) {
		date += earth_period;
	}
	return date;
}

function calculate_seasons() {
	let terms = [-pi/4, 0, pi/4, pi/2, 3*pi/4, pi, -3*pi/4, -pi/2, -pi/4, 0]
	let dates = [];
	let astro_seasons = [];
	let solar_seasons = [];
	for (let i = 0; i < terms.length; i++) {
		dates.push(date_from_longitude(terms[i]));
	}
	for (let i = 0; i < terms.length - 2; i += 2) {
		var astro_season_length = dates[i + 3] - dates[i + 1];
		if (astro_season_length < 0) {
			astro_season_length += earth_period;
		}
		var solar_season_length = dates[i + 2] - dates[i];
		if (solar_season_length < 0) {
			solar_season_length += earth_period;
		}
		astro_seasons.push(astro_season_length);
		solar_seasons.push(solar_season_length);

	}
	return [dates, astro_seasons, solar_seasons];
}

function calculate_precessions(m) {
	let nodal = -3/4*m + 9/32*Math.pow(m,2) + 273/128*Math.pow(m,3) + 9797/2048*Math.pow(m,4) + 199273/24576*Math.pow(m,5) + 6657733/589825*Math.pow(m,6);
	let apsidal = 3/4*m + 225/32*Math.pow(m,2) + 4071/128*Math.pow(m,3) + 265493/2048*Math.pow(m,4) + 12822631/24576*Math.pow(m,5) + 1273925965/589824*Math.pow(m,6) + 66702631253/7077888*Math.pow(m,7) + 29726828924189/679477248*Math.pow(m,8);
	nodal = 	1/nodal * earth_period;
	apsidal = 1/apsidal * earth_period;

	nodal = -6793.47709618;
	apsidal = 3233.0;

	let draconic = -nodal * moon_period / (-nodal + moon_period);
	let anomalistic = apsidal * moon_period / (apsidal - moon_period);
	let eclipse_year = -nodal * earth_period / (-nodal + earth_period);
	return [nodal, apsidal , draconic, anomalistic, eclipse_year];
}
