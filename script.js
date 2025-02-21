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
<tr  class="title"> <td colspan="3"> Orbital Properties </td> </tr>
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
<tr class="even">
	<td> Apoapsis (km) </td>
	<td> ${earth_apoapsis} </td>
	<td> ${moon_apoapsis} </td>
</tr>
<tr>
	<td> Sidereal Period (days) </td>
	<td> ${earth_period} </td>
	<td> ${moon_period} </td>
</tr>
<tr class="even">
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
let text = `
<thead>
<tr  class="title"> <td colspan="2"> Solar Calendar Properties </td> </tr>
<tr> 
	<td> Days in Year </td>
	<td> ${Math.floor(earth_period)} </td>
</tr>
<tr > <td colspan="2"> </td> </tr>
<tr> 
	<td> Leap Days per Years </td>
	<td> Error (years/1day) </td>
</tr>
</thead>
<tbody>
`;

var solar_leaps = dec_to_frac(earth_period, 10);
for (let i = 0; i < solar_leaps.length; i++) {
	text += "<tr ";
	if (i % 2 == 1) {
		text += "class='even'"
	}
	text += "> <td>" + solar_leaps[i][0].toString() + " / " + solar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(earth_period) + solar_leaps[i][0]/solar_leaps[i][1] - earth_period)).toString() + "</td> </tr>";
	solar_calendar_table.innerHTML = text;
}
solar_calendar_table.innerHTML += "</tbody>";
body.appendChild(solar_calendar_table);

const lunar_calendar_table = document.createElement("table");
text = `
<thead>
<tr  class="title"> <td colspan="2"> Lunar Calendar Properties </td> </tr>
<tr > 
	<td> Days in Short Month </td>
	<td> ${Math.floor(moon_synodic)} </td>
</tr>
<tr > 
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

var lunar_leaps = dec_to_frac(moon_synodic, 10);
for (let i = 0; i < lunar_leaps.length; i++) {
	text += "<tr ";
	if (i % 2 == 1) {
		text += "class='even'"
	}
	text += "> <td>" + lunar_leaps[i][0].toString() + " / " + lunar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(moon_synodic) + lunar_leaps[i][0]/lunar_leaps[i][1] - moon_synodic)).toString() + "</td> </tr>";
	lunar_calendar_table.innerHTML = text;
}

lunar_calendar_table.innerHTML += "</tbody>";
body.appendChild(lunar_calendar_table);

const lunisolar_calendar_table = document.createElement("table");
text = `
<thead>
<tr  class="title"> <td colspan="2"> Lunisolar Calendar Properties </td> </tr>
<tr > 
	<td> Days in Short Month </td>
	<td> ${Math.floor(moon_synodic)} </td>
</tr>
<tr > 
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

var lunisolar_leaps = dec_to_frac(months_in_year, 10);
for (let i = 0; i < lunisolar_leaps.length; i++) {
	text += "<tr ";
	if (i % 2 == 1) {
		text += "class='even'"
	}
	text += "> <td>" + lunisolar_leaps[i][0].toString() + " / " + lunisolar_leaps[i][1].toString() + "</td> <td>" + (1/(Math.floor(months_in_year) + lunisolar_leaps[i][0]/lunisolar_leaps[i][1] - months_in_year)).toString() + "</td> </tr>";
	lunisolar_calendar_table.innerHTML = text;
}
lunisolar_calendar_table.innerHTML += "</tbody>";
body.appendChild(lunisolar_calendar_table);

body.innerHTML += "<h2> Seasonal Properties </h2>";
var seasons = calculate_seasons();

const solar_terms_table = document.createElement("table");
solar_terms_table.innerHTML = `
<thead>
<tr  class="title"> <td colspan="3"> Solar Terms </td> </tr>
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
<tr class='even'> 
	<td> Spring Equinox </td>
	<td> ${seasons[0][1]} </td>
	<td> Autumn Equinox </td>
</tr>
<tr> 
	<td> Lixia / Bealtaine </td>
	<td> ${seasons[0][2]} </td>
	<td> Lidong / Samhain </td>
</tr>
<tr class='even'> 
	<td> Summer Solstice </td>
	<td> ${seasons[0][3]} </td>
	<td> Winter Solstice </td>
</tr>
<tr> 
	<td> Liqiu / Lunasa </td>
	<td> ${seasons[0][4]} </td>
	<td> Lichun / Imbolc </td>
</tr>
<tr class='even'> 
	<td> Autumn Equinox </td>
	<td> ${seasons[0][5]} </td>
	<td> Spring Equinox </td>
</tr>
<tr> 
	<td> Lidong / Samhain </td>
	<td> ${seasons[0][6]} </td>
	<td> Lixia / Bealtaine </td>
</tr>
<tr class='even'> 
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
<tr  class="title"> <td colspan="3"> Lengths of Astronomical Seasons </td> </tr>
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
<tr class='even'> 
	<td> Summer  </td>
	<td> ${seasons[1][1]} </td>
	<td> Winter </td>
</tr>
<tr> 
	<td> Autumn </td>
	<td> ${seasons[1][2]} </td>
	<td> Spring </td>
</tr>
<tr class='even'> 
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
<tr  class="title"> <td colspan="3"> Lengths of Solar Seasons </td> </tr>
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
<tr class='even'> 
	<td> Summer  </td>
	<td> ${seasons[2][1]} </td>
	<td> Winter </td>
</tr>
<tr> 
	<td> Autumn </td>
	<td> ${seasons[2][2]} </td>
	<td> Spring </td>
</tr>
<tr class='even'> 
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
<tr  class="title"> <td colspan="3"> Precession Properties </td> </tr>
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
<tr class='even'> 
	<td> Lunar Apsidal Precession </td>
	<td> ${precessions[1]} </td>
	<td> days / rev </td>
</tr>
<tr> 
	<td> Anomalistic Month </td>
	<td> ${precessions[2]} </td>
	<td> days </td>
</tr>
<tr class='even'> 
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
text = `
<thead>
<tr  class="title"> <td colspan="5"> Eclipse Cycles </td> </tr>
<tr> 
	<td> SM : half-DM </td>
	<td> Delta xi (deg) </td>
	<td> Planet Years </td>
	<td> Anom. Months </td>
	<td> Node </td>
</tr>

</thead>
<tbody>
`;

var eclipse_cycles = dec_to_frac((precessions[3]/2) / moon_synodic , 12);
for (let i = 0; i < eclipse_cycles.length; i++) {
	let quotient = eclipse_cycles[i][0] * moon_synodic / precessions[4];
	let error = (quotient - Math.floor(quotient)) * 360;
	let error_2 = (quotient - Math.floor(quotient) - 1) * 360;
	let error_3 = (quotient - Math.floor(quotient) - 0.5) * 360;
	if (Math.abs(error_2) < Math.abs(error)) {
		error = error_2;
	}
	if (Math.abs(error_3) < Math.abs(error)) {
		error = error_3;
	}
	text += "<tr "
	if (i % 2 == 1) {
		text += "class='even'"
	}
	text += "> <td>" + eclipse_cycles[i][0].toString() + " : " + eclipse_cycles[i][1].toString() + "</td> <td>" + error.toString() + "</td> ";
	text += "<td>" + (eclipse_cycles[i][0] * moon_synodic / earth_period).toString() + "</td> <td>" + (eclipse_cycles[i][0] * moon_synodic / precessions[2]).toString() + "</td> ";
	if (eclipse_cycles[i][1] % 2 == 0) {
		text += "<td> Same </td> </tr>"
	}
	else {
		text += "<td> Alternate </td> </tr>"
	}
	eclipse_cycle_table.innerHTML = text;
}

eclipse_cycle_table.innerHTML += "</tbody>";
body.appendChild(eclipse_cycle_table);

const eclipse_seasons_table = document.createElement("table");
eclipse_seasons = calculate_eclipse_seasons();
eclipse_seasons_table.innerHTML = `
<thead>
<tr  class="title"> <td colspan="4"> Average Eclipse Conditions </td> </tr>
<tr> 
	<td> Eclipse Type </td>
	<td> Beta less than (deg) </td>
	<td> Xi less than (deg) </td>
	<td> Season Length (days) </td>
</tr>
</thead>
<tbody>
<tr> 
	<td> Partial Solar </td>
	<td> ${eclipse_seasons[0][0] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][0] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][0] / (pi / earth_period)} </td>
</tr>
<tr class='even'> 
	<td> Central Solar </td>
	<td> ${eclipse_seasons[0][1] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][1] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][1] / (pi / earth_period)} </td>
</tr>
<tr> 
	<td> Penumbral Lunar </td>
	<td> ${eclipse_seasons[0][2] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][2] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][2] / (pi / earth_period)} </td>
</tr>
<tr class='even'> 
	<td> Partial Lunar </td>
	<td> ${eclipse_seasons[0][3] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][3] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][3] / (pi / earth_period)} </td>
</tr>
<tr> 
	<td> Total Lunar </td>
	<td> ${eclipse_seasons[0][4] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][4] * 180 / pi} </td>
	<td> ${eclipse_seasons[1][4] / (pi / earth_period)} </td>
</tr>
</tbody>
`;
body.appendChild(eclipse_seasons_table);

const eclipse_cycle_longevity_table = document.createElement("table");
text = `
<thead>
<tr  class="title"> <td colspan="6"> Approximate Eclipse Cycle Longevities </td> </tr>
<tr> 
	<td> Eclipse Cycle </td>
	<td> Total Num Solar </td>
	<td> Num Central Solar </td>
	<td> Total Num Lunar </td>
	<td> Num Total Lunar </td>
	<td> Longevity (years) </td>
</tr>

</thead>
<tbody>
`;
for (let i = 0; i < eclipse_cycles.length; i++) {
	let quotient = eclipse_cycles[i][0] * moon_synodic / precessions[4];
	let error = (quotient - Math.floor(quotient)) * 2 * pi;
	let error_2 = (quotient - Math.floor(quotient) - 1) * 2 * pi;
	let error_3 = (quotient - Math.floor(quotient) - 0.5) * 2 * pi;
	if (Math.abs(error_2) < Math.abs(error)) {
		error = error_2;
	}
	if (Math.abs(error_3) < Math.abs(error)) {
		error = error_3;
	}
	error = Math.abs(error);
	text += "<tr "
	if (i % 2 == 1) {
		text += "class='even'"
	}
	text += "> <td>" + eclipse_cycles[i][0].toString() + " : " + eclipse_cycles[i][1].toString() + "</td>"
	for (let j = 0; j < eclipse_seasons[1].length; j++) {
		if (j != 3) {
			text +=  "<td>" + (2 + Math.floor(2 * eclipse_seasons[1][j] / error)).toString() + "</td>";
		}
	}
	text += "<td>" + (eclipse_cycles[i][0] * moon_synodic / earth_period) * (1 + Math.floor(2 * eclipse_seasons[1][0] / error)).toString() + " </td></tr>";

	eclipse_cycle_longevity_table.innerHTML = text;
}

eclipse_cycle_longevity_table.innerHTML += "</tbody>";
body.appendChild(eclipse_cycle_longevity_table);

function period(mass1, mass2, sma) {
	return Math.sqrt(4 * pi * pi * sma * sma * sma / (G * (mass1 + mass2))) / 3600 / day;
}

function surface_gravity(M, R) {
	return G * M / R / R * 1000;
}

function dec_to_frac(number, n) {
	let num_approxs = n;
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

function calculate_eclipse_seasons() {
	let moon_semidiameter = Math.asin(moon_radius / moon_sma);
	let sun_semidiameter = Math.asin(sun_radius / earth_sma);
	let moon_parallax = Math.asin(earth_radius / moon_sma);
	let sun_parallax = Math.asin(earth_radius / earth_sma);
	let moon_movement = 2 * pi / moon_period;
	let sun_movement = 2 * pi / earth_period;
	let q = moon_movement / sun_movement;
	let i_prime = Math.atan(q / (q - 1) * Math.tan(moon_inclination));

	let partial_solar = 1/Math.cos(i_prime) * (moon_semidiameter + sun_semidiameter + moon_parallax - sun_parallax);
	let central_solar = 1/Math.cos(i_prime) * (moon_semidiameter - sun_semidiameter + moon_parallax - sun_parallax);
	let penumbral_lunar = 1/Math.cos(i_prime) * (moon_semidiameter + sun_semidiameter + moon_parallax + sun_parallax);
	let partial_lunar = 1/Math.cos(i_prime) * (moon_semidiameter - sun_semidiameter + moon_parallax - sun_parallax);
	let total_lunar = 1/Math.cos(i_prime) * (-moon_semidiameter - sun_semidiameter + moon_parallax - sun_parallax);
	
	let beta_array = [partial_solar, central_solar, penumbral_lunar, partial_lunar, total_lunar];
	let xi_array = [];
	for (let i = 0; i < beta_array.length; i++) {
		xi_array.push(Math.asin(Math.tan(beta_array[i]) / Math.tan(moon_inclination)));
	}
	return [beta_array, xi_array];
}
