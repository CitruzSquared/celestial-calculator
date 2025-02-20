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

body.innerHTML += "<h2> Basic Properties </h2>";
const basic_table = document.createElement("table");
basic_table.innerHTML = `
<thead>
<tr> <td colspan="3"> Orbital Properties </td> </tr>
<tr> 
	<td> Property </td>
	<td> Earth </td>
	<td> Moon </td>
</tr>
</thead>
<tbody>
<tr>
	<td> Periapsis (km) </td>
	<td> ${earth_periapsis} </td>
	<td> ${moon_periapsis} </td>
</tr>
<tr>
	<td> Apoapsis (km) </td>
	<td> ${earth_apoapsis} </td>
	<td> ${moon_apoapsis} </td>
</tr>
<tr>
	<td> Sidereal Period (days) </td>
	<td> ${earth_period} </td>
	<td> ${moon_period} </td>
</tr>
<tr>
	<td> Synodic Period (days) </td>
	<td> - </td>
	<td> ${moon_synodic} </td>
</tr>
<tr>
	<td> Surface Gravity (m/s/s) </td>
	<td> ${surface_gravity(earth_mass, earth_radius)} </td>
	<td> ${surface_gravity(moon_mass, moon_radius)} </td>
</tr>
</tbody>
`;
body.appendChild(basic_table);

body.innerHTML += "<h2> Calendar Properties </h2>";
const solar_calendar_table = document.createElement("table");
solar_calendar_table.innerHTML = `
<thead>
<tr> <td colspan="2"> Solar Calendar Properties </td> </tr>
<tr> 
	<td> Days in Year </td>
	<td> ${Math.floor(earth_period)} </td>
</tr>
<tr> <td colspan="2"> </td> </tr>
<tr> 
	<td> Leap Days per Years </td>
	<td> Error (years/1day) </td>
</tr>
</thead>
<tbody>
`;

var solar_leaps = dec_to_frac(earth_period);
for (let i = 0; i < solar_leaps.length; i++) {
	solar_calendar_table.innerHTML += "<tr> <td>" + solar_leaps[i][0].toString() + " / " + solar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(earth_period) + solar_leaps[i][0]/solar_leaps[i][1] - earth_period)).toString() + "</td> </tr>";
}
solar_calendar_table.innerHTML += "</tbody>";
body.appendChild(solar_calendar_table);

const lunar_calendar_table = document.createElement("table");
lunar_calendar_table.innerHTML = `
<thead>
<tr> <td colspan="2"> Lunar Calendar Properties </td> </tr>
<tr> 
	<td> Days in Short Month </td>
	<td> ${Math.floor(moon_synodic)} </td>
</tr>
<tr> 
	<td> Months in Year </td>
	<td> ${Math.floor(months_in_year)} </td>
</tr>
<tr> 
	<td> Drift per Year (days) </td>
	<td> ${moon_synodic * Math.floor(months_in_year) - earth_period} </td>
</tr>
<tr> <td colspan="2"> </td> </tr>
<tr> 
	<td> Long Months per Months </td>
	<td> Error (months/1day) </td>
</tr>
</thead>
<tbody>
`;

var lunar_leaps = dec_to_frac(moon_synodic);
for (let i = 0; i < lunar_leaps.length; i++) {
	lunar_calendar_table.innerHTML += "<tr> <td>" + lunar_leaps[i][0].toString() + " / " + lunar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(moon_synodic) + lunar_leaps[i][0]/lunar_leaps[i][1] - moon_synodic)).toString() + "</td> </tr>";
}

lunar_calendar_table.innerHTML += "</tbody>";
body.appendChild(lunar_calendar_table);

const lunisolar_calendar_table = document.createElement("table");
lunisolar_calendar_table.innerHTML = `
<thead>
<tr> <td colspan="2"> Lunisolar Calendar Properties </td> </tr>
<tr> 
	<td> Days in Short Month </td>
	<td> ${Math.floor(moon_synodic)} </td>
</tr>
<tr> 
	<td> Months in Year </td>
	<td> ${Math.floor(months_in_year)} </td>
</tr>
<tr> <td colspan="2"> </td> </tr>
<tr> 
	<td> Leap Months per Years </td>
	<td> Error (years/1month) </td>
</tr>
</thead>
<tbody>
`;

var lunisolar_leaps = dec_to_frac(months_in_year);
for (let i = 0; i < lunisolar_leaps.length; i++) {
	lunisolar_calendar_table.innerHTML += "<tr> <td>" + lunisolar_leaps[i][0].toString() + " / " + lunisolar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(months_in_year) + lunisolar_leaps[i][0]/lunisolar_leaps[i][1] - months_in_year)).toString() + "</td> </tr>";
}
lunisolar_calendar_table.innerHTML += "</tbody>";
body.appendChild(lunisolar_calendar_table);

body.innerHTML += "<h2> Seasonal Properties </h2>";
var seasons = calculate_seasons();

const solar_terms_table = document.createElement("table");
solar_terms_table.innerHTML = `
<thead>
<tr> <td colspan="3"> Solar Terms </td> </tr>
<tr> 
	<td> Northern Term </td>
	<td> Days after Periapsis </td>
	<td> Southern Term </td>
</tr>
</thead>
<tbody>
<tr> 
	<td> Lichun / Imbolc </td>
	<td> ${seasons[0][0]} </td>
	<td> Liqiu / Lunasa </td>
</tr>
<tr> 
	<td> Spring Equinox </td>
	<td> ${seasons[0][1]} </td>
	<td> Autumn Equinox </td>
</tr>
<tr> 
	<td> Lixia / Bealtaine </td>
	<td> ${seasons[0][2]} </td>
	<td> Lidong / Samhain </td>
</tr>
<tr> 
	<td> Summer Solstice </td>
	<td> ${seasons[0][3]} </td>
	<td> Winter Solstice </td>
</tr>
<tr> 
	<td> Liqiu / Lunasa </td>
	<td> ${seasons[0][4]} </td>
	<td> Lichun / Imbolc </td>
</tr>
<tr> 
	<td> Autumn Equinox </td>
	<td> ${seasons[0][5]} </td>
	<td> Spring Equinox </td>
</tr>
<tr> 
	<td> Lidong / Samhain </td>
	<td> ${seasons[0][6]} </td>
	<td> Lixia / Bealtaine </td>
</tr>
<tr> 
	<td> Winter Solstice </td>
	<td> ${seasons[0][7]} </td>
	<td> Summer Solstice </td>
</tr>
</tbody>
`;
body.appendChild(solar_terms_table);

const astro_seasons_table = document.createElement("table");
astro_seasons_table.innerHTML = `
<thead>
<tr> <td colspan="3"> Lengths of Astronomical Seasons </td> </tr>
<tr> 
	<td> Northern Season </td>
	<td> Length (days) </td>
	<td> Southern Season </td>
</tr>
</thead>
<tbody>
<tr> 
	<td> Spring </td>
	<td> ${seasons[1][0]} </td>
	<td> Autumn </td>
</tr>
<tr> 
	<td> Summer  </td>
	<td> ${seasons[1][1]} </td>
	<td> Winter </td>
</tr>
<tr> 
	<td> Autumn </td>
	<td> ${seasons[1][2]} </td>
	<td> Spring </td>
</tr>
<tr> 
	<td> Winter </td>
	<td> ${seasons[1][3]} </td>
	<td> Summer  </td>
</tr>
</tbody>
`;
body.appendChild(astro_seasons_table);


const solar_seasons_table = document.createElement("table");
solar_seasons_table.innerHTML = `
<thead>
<tr> <td colspan="3"> Lengths of Solar Seasons </td> </tr>
<tr> 
	<td> Northern Season </td>
	<td> Length (days) </td>
	<td> Southern Season </td>
</tr>
</thead>
<tbody>
<tr> 
	<td> Spring </td>
	<td> ${seasons[2][0]} </td>
	<td> Autumn </td>
</tr>
<tr> 
	<td> Summer  </td>
	<td> ${seasons[2][1]} </td>
	<td> Winter </td>
</tr>
<tr> 
	<td> Autumn </td>
	<td> ${seasons[2][2]} </td>
	<td> Spring </td>
</tr>
<tr> 
	<td> Winter </td>
	<td> ${seasons[2][3]} </td>
	<td> Summer  </td>
</tr>
</tbody>
`;
body.appendChild(solar_seasons_table);

body.innerHTML += "<h2> Eclipse Properties </h2>";
var precessions = calculate_precessions(m);

const precession_table = document.createElement("table");
precession_table.innerHTML = `
<thead>
<tr> <td colspan="3"> Precession Properties </td> </tr>
<tr> 
	<td> Property </td>
	<td> Value </td>
	<td> Unit </td>
</tr>
</thead>
<tbody>
<tr> 
	<td> Lunar Nodal Precession </td>
	<td> ${precessions[0]} </td>
	<td> days / rev </td>
</tr>
<tr> 
	<td> Lunar Apsidal Precession </td>
	<td> ${precessions[1]} </td>
	<td> days / rev </td>
</tr>
<tr> 
	<td> Anomalistic Month </td>
	<td> ${precessions[2]} </td>
	<td> days </td>
</tr>
<tr> 
	<td> Draconic Month </td>
	<td> ${precessions[3]} </td>
	<td> days </td>
</tr>
<tr> 
	<td> Eclipse Year </td>
	<td> ${precessions[4]} </td>
	<td> days </td>
</tr>
</tbody>
`;
body.appendChild(precession_table);

const eclipse_cycle_table = document.createElement("table");
eclipse_cycle_table.innerHTML = `
<thead>
<tr> <td colspan="4"> Eclipse Cycles </td> </tr>
<tr> 
	<td> Half-EY : SM </td>
	<td> Delta xi (degrees) </td>
	<td> Years </td>
	<td> A.M. </td>
</tr>

</thead>
<tbody>
`;

var eclipse_cycles = dec_to_frac(moon_synodic/(precessions[4]/2));
for (let i = 0; i < eclipse_cycles.length; i++) {
	let quotient = eclipse_cycles[i][1] * moon_synodic / precessions[4];
	let error = (quotient - Math.floor(quotient)) * 360;
	error = Math.min(error, 360 - error);
	error = Math.min(error, Math.abs(error-180), Math.abs(180-error));
	let text = "<tr> <td>" + eclipse_cycles[i][0].toString() + " : " + eclipse_cycles[i][1].toString() + "</td> <td>" + error.toString() + "</td> ";
	text += "<td>" + (eclipse_cycles[i][1] * moon_synodic / earth_period).toString() + "</td> <td>" + (eclipse_cycles[i][1] * moon_synodic / precessions[2]).toString() + "</td> </tr>";
	eclipse_cycle_table.innerHTML += text;
}

eclipse_cycle_table.innerHTML += "</tbody>";
body.appendChild(eclipse_cycle_table);

body.innerHTML += "<p> Eclipse cycles: " + dec_to_frac(moon_synodic/(precessions[4]/2)).join(" or ") + " half-E.Y./S.M.</p>";

function period(mass1, mass2, sma) {
	return Math.sqrt(4 * pi * pi * sma * sma * sma / (G * (mass1 + mass2))) / 3600 / day;
}

function surface_gravity(M, R) {
	return G * M / R / R * 1000;
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
		
		fracs.push([numerator,denominator]);
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
	apsidal = 3232.6054285

	let anomalistic = apsidal * moon_period / (apsidal - moon_period);
	let draconic = -nodal * moon_period / (-nodal + moon_period);
	let eclipse_year = -nodal * earth_period / (-nodal + earth_period);
	return [nodal, apsidal , anomalistic, draconic, eclipse_year];
}
