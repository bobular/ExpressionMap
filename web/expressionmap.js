/*
 * requires: prototype, raphael, scriptaculous
 *
 */

//PUT THESE INTO config.js
var fold_thresh = 1.5;
var num_conds = 5;
var z = 1000;  // a global to keep popups always on top
var nodepapers = new Array(); // a global to keep svg 'paper' refs in (NOTE ORDER: [y][x])
var tempelements = new Array(); // [class]=Array keep track of svg elements which may need to be removed
var sliders = new Array();
var params = document.location.search.toQueryParams();

function em_draw_nodes(map_id, map_width, map_height, maxgenes) {
	var sqrtmaxgenes = Math.sqrt(maxgenes);
	var sumgenes = 0;
	for (var y=0; y<map_height; y++) {
		nodepapers[y] = new Array();
    for (var x=0; x<map_width; x++) {
			// add handlers to table elements
			var td = em_get_node(x,y);
			// td.update();
			var numgenes = td.readAttribute("numgenes");
			sumgenes += parseInt(numgenes);
			var paper = Raphael(td, 31, 31);
			paper.circle(15, 15, 15*Math.sqrt(numgenes)/sqrtmaxgenes).
				attr({ "fill":"#ccc", "stroke-width":0, "title": "node "+x+","+y+" ("+numgenes+" genes) click for more..." });
			nodepapers[y][x] = paper; // for future reference (hope it's not too memory intensive)
		}
	}
	$('please-wait').hide();
	$('num-genes').update("Genes: <b>"+sumgenes+"</b>");
}


function em_activate_nodes(map_id, map_width, map_height) {
	for (var y=0; y<map_height; y++) {
    for (var x=0; x<map_width; x++) {
			// add handlers to table elements
			var td = em_get_node(x,y);
			var numgenes = td.readAttribute("numgenes");
			var popup = new Element('div', { id: 'popup-'+x+'-'+y }).addClassName('node-popup').hide();
			td.insert({ top: popup });
			var header = new Element('div').addClassName('node-popup-header').update("&nbsp; node "+x+","+y+" ("+numgenes+" genes)");
			popup.insert({ top: header });
			var xclose = new Element('img', { src: 'close.png', width: 15, height: 15, align: 'left' });
			xclose.observe('click', em_node_hide(x,y));
			header.insert({ top: xclose });

      // try to add something to the paper
			/*
			var paper = nodepapers[y][x];
			var circle = paper.circle(15, 15, 2).
				attr({ "fill":"red", "stroke-width":0 });
			tempelements.reddots.push(circle);
			*/

			td.observe('click', em_node_callback(x,y));
		}

		
	}

	// removing tempelements by key
	/*
	setTimeout(function() {
			tempelements.reddots.each(function(element) {
					element.remove();
				});
		}, 3000);
	*/

  function em_node_callback(x,y) {
		return function(event) {
			em_get_node(x,y).stopObserving();
			em_display_node(map_id, x, y);
		};
	}
  function em_node_hide(x,y) {
		return function(event) {
			var td = em_get_node(x,y);
			var popup = td.down('.node-popup');
			// only works with defer or a timeout/delay
			Element.hide.defer(popup);
			// ditto...
			setTimeout(function(){ td.observe('click', em_node_callback(x,y)); }, 100);
		};
	}

}

function em_display_node(map_id, x, y) {
	var td = em_get_node(x,y);
	var popup = td.down('div.node-popup');
	popup.setStyle({zIndex: z++ });
	popup.show();
	new Draggable(popup.id, { handle: 'node-popup-header' });
	var info = popup.down('.info');
	if (info) {
		// already there so do nothing
	} else {
		info = new Element('div').addClassName('info');

		info.insert({ bottom: new Element('div').addClassName("loading-data").update("<br/>loading data...") });
		popup.insert({bottom: info});

		// by threshold:
		// new Ajax.Request('json.cgi?action=nodeinfo&x='+x+'&y='+y+'&id='+map_id+'&ft='+fold_thresh, {
		// by num_conds most changed:
		new Ajax.Request('json.cgi?action=nodeinfo&x='+x+'&y='+y+'&id='+map_id+'&nc='+num_conds, {
				method:'get',
					onSuccess: function(response){ em_handle_node_info(info, response.responseText.evalJSON()) },
					onFailure: function(response){ info.update('Sorry, an error occurred - please email the help desk.') }
				});
	}

	// clear any search result divs
  info.select("div.search-results").each(function(div) { div.remove() });

	// show the search results if there are any
	if (td.searchResults && td.searchResults.size()) {
		for (var search_num=3; search_num>=1; search_num--) {
			if (td.searchResults[search_num]) {
				var result_div = new Element('div').addClassName('search-results');
				result_div.insert({ top: new Element('div').addClassName('popup-header').addClassName('search'+search_num).update("Query "+search_num+":") });
				result_div.insert({ bottom: new Element('div').addClassName('search'+search_num).update(td.searchResults[search_num].collect(em_gene_link).join(", ")) });
				info.insert({ top: result_div });
			}
		}
		info.insert({ top: new Element('div').addClassName('node-popup-subheader').addClassName('search-results').update("Gene query results:") });
	}


}

function em_handle_node_info(element, data) {
	element.down(".loading-data").remove();
	element.insert({bottom: new Element('div').addClassName('node-popup-subheader').update('Expression characteristics ('+num_conds+' largest changes):')});
	var weights = new Element('div').addClassName('weights');
	var keys = Object.keys(data.weights);
	if (keys.size()) {
		var table = new Element('table');
		var last_weight;
		keys.sortBy(function(key){return -parseFloat(data.weights[key]);}).each(function(key){
				var weight = data.weights[key];
				var row = new Element('tr');
				var hue = 120-weight*50;
				if (hue>240) { hue=240 }
				if (hue<0) { hue = 0 }
				var col = Raphael.hsl(hue, 100, 40);
				row.setStyle({ color: col });
				row.insert({bottom: new Element('td').update(weight > 0 ? Math.pow(2,weight).toFixed(1)+'-fold&nbsp;up' : Math.pow(2,-weight).toFixed(1)+'-fold&nbsp;down')});
				row.insert({bottom: new Element('td').update(key)});
				if (last_weight !== null && last_weight > 0 && weight < 0) {
					row.addClassName('top-rule');
				}
				table.insert({bottom: row});
				last_weight = weight;
			});
		weights.insert({bottom: table});
	} else {
		weights.update('none');
	}
	element.insert({bottom: weights});

	element.insert({bottom: new Element('div').addClassName('node-popup-subheader').update('Genes:')});
	element.insert({bottom: new Element('div').
				addClassName('genes').
				update(data.genes.collect(em_gene_link).join(", "))});
	element.insert({bottom: new Element('div').addClassName('node-popup-subheader').update('Over-represented Gene Ontology terms:')});
	var goterm_table = new Element('table');
	
	element.insert({bottom: new Element('div').addClassName('goterms').insert({bottom: goterm_table })});
		
	if (data.go_ora_results.size()) {
		var goterm_headings = new Element('tr');
		goterm_table.insert({bottom: goterm_headings});
		goterm_headings.insert({bottom: new Element('th').update('#genes') });
		goterm_headings.insert({bottom: new Element('th').update('Accession') });
		goterm_headings.insert({bottom: new Element('th').update('Name') });
		goterm_headings.insert({bottom: new Element('th').update('Type') });
		goterm_headings.insert({bottom: new Element('th').update('p-value') });

		data.go_ora_results.each(function(result){
				var goterm_row = new Element('tr');
				goterm_table.insert({bottom: goterm_row});
				goterm_row.insert({bottom: new Element('td').update(result.score) });
				goterm_row.insert({bottom: new Element('td').update(result.acc) });
				goterm_row.insert({bottom: new Element('td').update(result.term) });
				goterm_row.insert({bottom: new Element('td').update(result.type) });
				goterm_row.insert({bottom: new Element('td').update(result.pvalue) });
			});
	} else {
		var goterm_row = new Element('tr');
		goterm_table.insert({bottom: goterm_row});
		goterm_row.insert({bottom: new Element('td').update("none with corrected p-value < 0.05") });
	}
}


/*
 * returns the <td> given x,y
 */

function em_get_node(x, y) {
  return $('node-'+x+'-'+y);
}


/*
 * em_ehighlight - run the expression filter using form data
 */

function em_ehighlight(map_id, filter_num) {
	// clear any previous highlighting for this filter
	if (tempelements['ehighlight'+filter_num] == null) {
		tempelements['ehighlight'+filter_num] = new Array();
	}
	var temp = tempelements['ehighlight'+filter_num];
	temp.each(function(element) {
					element.remove();
		});

	em_hide_popups();

	// send the ajax request
	var form = $('control-panel');
	var vecdim_idx = form['e'+filter_num].getValue();
  var log2_thresh = $('efilterval'+filter_num).sliderValue; // Form.getInputs('control-panel','radio','t'+filter_num).find(function(radio) { return radio.checked; }).value;
	if (vecdim_idx > -1 && log2_thresh != 0) {
		var effect = Effect.Pulsate('map1', { duration: 5.0, pulses: 2 });
		new Ajax.Request('json.cgi?action=ehighlight&v='+vecdim_idx+'&lt='+log2_thresh+'&id='+map_id, {
				method:'get',
					onSuccess: function(response){ em_handle_highlight_response(filter_num, response.responseText.evalJSON(), effect) },
					onFailure: function(response){ alert('Sorry, an error occurred - please email the help desk.') }
		});
	}
}

function em_handle_highlight_response(filter_num, data, effect) {
	var temp = tempelements['ehighlight'+filter_num];
	effect.cancel();
	Effect.Appear('map1', { duration: 0.1 });
	var nodes = data.nodes;
	nodes.each(function(node_coords) {
			var x = node_coords[0];
			var y = node_coords[1];
			var paper = nodepapers[y][x];
			if (filter_num == 1) {
				highlight1(paper, temp);
			} else if (filter_num == 2) {
				highlight2(paper, temp);
			} else if (filter_num == 3) {
				highlight3(paper, temp);
			}
		});
	$('efilternodes'+filter_num).innerHTML = ' &rarr; '+nodes.size()+' nodes';
}


function highlight1(paper, temp) {
	var path = paper.path("M 2 0 L 2 31 M 28 0 L 28 31 M 0 2 L 31 2 M 0 28 L 31 28").
		attr({ "stroke":"red", "stroke-width":2 });
	temp.push(path);
}
function highlight2(paper, temp) {
	var path = paper.path("M 4 0 L 4 31 M 26 0 L 26 31 M 0 4 L 31 4 M 0 26 L 31 26").
		attr({ "stroke":"green", "stroke-width":2 });
	temp.push(path);
}
function highlight3(paper, temp) {
	var path = paper.path("M 6 0 L 6 31 M 24 0 L 24 31 M 0 6 L 31 6 M 0 24 L 31 24").
		attr({ "stroke":"blue", "stroke-width":2 });
	temp.push(path);
}


function em_init_sliders(map_id) {
  for (var i=1; i<=3; i++) {
		var threshold = params['t'+i] || 0;
		var slide_handler = em_handle_slider_slide(i);
		var slider = new Control.Slider('handle'+i , 'track'+i, {
				  range: $R(-4,4),
					values: $R(-40,40).collect(function(val){ return val/10; }),
					onChange: em_handle_slider_change(map_id,i),
					onSlide: slide_handler
			});
  	slider.setValue(threshold);
    sliders.push(slider);
 	}
}

function em_init_autocompleters(map_id) {
	new Ajax.Request('json.cgi?action=conditions&id='+map_id, {
			method:'get',
				onSuccess: function(response){ em_handle_autocomplete_response(map_id, response.responseText.evalJSON()) },
				onFailure: function(response){ alert('Sorry, an error occurred - please email the help desk.') }
	});
}

function em_handle_autocomplete_response(map_id, data) {
	var index2condition = $H(data.conditions);
	var condition2index = $H();
	var conditions = index2condition.values();
	$('num-conditions').update("Conditions: <b>"+conditions.size()+"</b>");

	index2condition.each( function(pair) {
			condition2index[pair.value] = pair.key;
		});
	for (var i=1; i<=3; i++) {
    $('efiltertext'+i).insert({ after:
				new Element('div', {id:"efilterauto"+i, class:'autocomplete'})
					});

	  new Autocompleter.Local('efiltertext'+i,
														'efilterauto'+i,
														conditions,
														{ afterUpdateElement: em_autocomplete_change_menu(map_id, i, condition2index),
															choices: 20,
														});
	}
}

function em_autocomplete_change_menu(map_id, filter_num, condition2index) {
	return function(inputfield,listitem) {
		var index = condition2index[listitem.textContent];
		$('e'+filter_num).setValue(index);
		em_ehighlight(map_id,filter_num);
		Effect.Pulsate('handle'+filter_num, { duration: 1.0 });
	};
}


function em_handle_slider_change(map_id, slider_num) {
	return function(v){
		$('efilterval'+slider_num).innerHTML = v < 0 ? '&lt;'+v : v > 0 ? '&gt;'+v : 'no highlighting';
		$('efilterval'+slider_num).sliderValue = v;
		$('efilternodes'+slider_num).update(v == 0 ? '' : '(log<sub>2</sub> scale)');
		em_ehighlight(map_id,slider_num);
  	em_update_permalink($('control-panel'));
	};
}

function em_handle_slider_slide(slider_num) {
	return function(v){
		$('efilterval'+slider_num).innerHTML =  v < 0 ? '&lt;'+v : v > 0 ? '&gt;'+v : 'no highlighting';
		$('efilterval'+slider_num).sliderValue = v;
		$('efilternodes'+slider_num).update(v == 0 ? '' : '(log<sub>2</sub> scale)');
	};
}

function em_handle_condition_change(map_id, menu, filter_num) {
  var val=menu.value;
	if (val > -1) {
		$('efiltertext'+filter_num).value = $A(menu.options).filter(function(o){return o.value==val;}).first().label;
		Effect.Pulsate('handle'+filter_num, { duration: 1.0 });
	} else {
		$('efiltertext'+filter_num).value = '';
		$('efilterval'+filter_num).update();
		$('efilternodes'+filter_num).update();
		sliders[filter_num].setValue(0);
	}
  em_ehighlight(map_id,filter_num);
}

function em_reset_all() {
	sliders.each(function(slider) { slider.setValue(0); });
	var form = $('control-panel');
  for (var filter_num=1; filter_num<=3; filter_num++) {
		$('efilterval'+filter_num).update('no highlighting');
		$('efilternodes'+filter_num).update();
		$('efiltertext'+filter_num).setValue("");
		$('e'+filter_num).setValue(-1);
	}
	for (var search_num=1; search_num<=3; search_num++) {
		var temp = tempelements['search'+search_num];
		if (temp) {
			temp.each(function(element) {
					element.remove();
				});
		}
		$('search'+search_num).setValue("");
	}
	$('map1').select('td').each(function(td) {
			if (td.searchResults) {
				delete td.searchResults;
			}
		});
	em_hide_popups();
	em_update_permalink(form);
}

function em_hide_popups() {
  $('map1').
		select('div.node-popup').
		filter(function(div) { return div.visible() }).
		each(function(div) {

				var img = div.down('img');

				try {
					img.click();
				} catch (e) {
					// doesn't work in Chrome
					// all hail http://blog.miranda.or.at/2011/07/javascript-why-simulating-a-click-using-the-click-method-doesnt-work-in-chrome !!
					var event = document.createEvent("MouseEvents");
					event.initMouseEvent("click", true, true, window,0, 0, 0, 0, 0, false, false, false, false, 0, null);
					var _r = !img.dispatchEvent(event);
				}
				//				var td = div.up('td');
				// td.observe('click', em_node_callback(td.cellIndex, td.up('tr').rowIndex));
			});
}

function em_handle_search_submit(map_id, search_num) {
	var query = $('search'+search_num).value;
	/* clear any previous labels for this search */
	if (tempelements['search'+search_num] == null) {
		tempelements['search'+search_num] = new Array();
	}
	var temp = tempelements['search'+search_num];
	temp.each(function(element) {
					element.remove();
		});

	$('map1').select('td').each(function(td) {
			if (td.searchResults && td.searchResults[search_num]) {
				delete td.searchResults[search_num];
			}
		});

	em_hide_popups();

	if (query.match(/.*\w.*/)) {
		var effect = Effect.Pulsate('map1', { duration: 5.0, pulses: 2 });
		new Ajax.Request('json.cgi?action=search&q='+escape(query)+'&id='+map_id, {
				method:'get',
					onSuccess: function(response){ em_handle_search_response(response.responseText.evalJSON(), search_num, effect) },
					onFailure: function(response){ alert('Sorry, an error occurred - please email the help desk.') }
		});
	}
}

function em_handle_search_response(data, search_num, effect) {
	/*	console.log(search_num);
			console.log(data);
	*/
	if (effect) {
		effect.cancel();
	}
	Effect.Appear('map1', { duration: 0.1 });


	/* now draw the search results */

	var temp = tempelements['search'+search_num];
	var nodes = data.search_results;
	nodes.each(function(node) {
			var x = node.x;
			var y = node.y;
			var paper = nodepapers[y][x];
			var colour;
			var lx = 15;
			var ly = 10;
			if (search_num == 1) {
				colour = 'yellow';
			} else if (search_num == 2) {
				colour = 'orange';
				ly = 15;
			} else if (search_num == 3) {
				colour = 'pink';
				ly = 20;
			}

			var hits = node.hits;
			var num = hits.size();
			var halfwidth = 14 * hits.size()/data.max_node_hits;
			if (halfwidth < 2.5) halfwidth = 2.5;
			var blob = paper.rect(lx-halfwidth, ly-2, 2*halfwidth, 5, 2).attr({ "stroke": "#333", "stroke-width": "1", "fill": colour, "title": 'Gene query '+search_num+' result: ' + hits.join(", ") + (num > 1 ? ' ('+num+' genes)' : '') });
			temp.push(blob);

			var td = em_get_node(x, y);
			if (td.searchResults == null) {
				td.searchResults = new Array();
			}
			td.searchResults[search_num] = hits;
		});
	
	if (nodes.size() == 0) {
		alert("No results for 'Gene query "+search_num+"'.  Please modify your query and try again.  TIPS: searches are whole-word only; no wildcards allowed; gene, GO and Interpro IDs must be in uppercase with the correct number of zeros.");
	}
}

function em_gene_link(string) {
	return string.replace(/^(\w+)/, '<a href="'+em_config.gene_expression_linkout+'$1">$1</a>');
}

function em_observe_form_for_permalink() {
	// do it once now
	em_update_permalink($('control-panel'));
	// and also any time the form changes
	new Form.Observer('control-panel', 0.67, function(form){
			em_update_permalink(form);
  });
}

function em_update_permalink(form) {
  var permalink = $('permalink');
	permalink.href = "?";
	var parvals = $A();
	if (form.search1.value) parvals.push("search1="+escape(form.search1.value));
	if (form.search2.value) parvals.push("search2="+escape(form.search2.value));
	if (form.search3.value) parvals.push("search3="+escape(form.search3.value));
	if (form.e1.value > -1) parvals.push("e1="+escape(form.e1.value), "t1="+($('efilterval1').sliderValue || 0));
	if (form.e2.value > -1) parvals.push("e2="+escape(form.e2.value), "t2="+($('efilterval2').sliderValue || 0));
	if (form.e3.value > -1) parvals.push("e3="+escape(form.e3.value), "t3="+($('efilterval3').sliderValue || 0));
	// join them together
	if (parvals.size()) { 
		permalink.href = "?"+parvals.join('&');
		permalink.show();
	} else {
		permalink.hide();
	}
}

function em_init_search_boxes(map_id) {
	// submit the searches if pre-filled with any data
	for (var search_num=1; search_num<=3; search_num++) {
		em_handle_search_submit(map_id, search_num);
	}
}
