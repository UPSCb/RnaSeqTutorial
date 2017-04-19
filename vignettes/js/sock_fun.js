/*
	We start by creating a hash out of some random things which will be the
	anonymous identifier for each participant. Most of the functions depend
	on the cookie `PARTID` being set correctly and having this hash as a value.
*/
	var PARTHASH = CryptoJS.MD5(Date() + Math.random()).toString();
	//Setting up socket.io
	var socket = io.connect('http://130.239.72.131:3001');
	//Quick function to set a cookie
	function setCookie(cname, cvalue, exdays) {
		var d = new Date();
		d.setTime(d.getTime() + (exdays*24*60*60*1000));
		var expires = "expires="+d.toUTCString();
		document.cookie = cname + "=" + cvalue + "; " + expires;
	}
	//The default lifetime of the participant cookie is 5 days.
	function setPartCookie() {
		setCookie("PARTID",PARTHASH, 5);
	}
	//We will also need to check if the cookie has already been set (to not overwrite)
	function checkPartCookie(){
	  if(document.cookie.match(/PARTID/g) != null){
		  return true;
	  } else {
		  return false;
	  }
	}
	//Function to extract the cookie
	function getPartCookie(){
		if(document.cookie.match(/PARTID=.+/g) != null){
		  var pcok = document.cookie.match(/PARTID=.+/g)[0].split(";")[0].replace("PARTID=","");
		  return pcok;
	    }
	}
	//Simle function to remove duplicate values from an array
	function uniq(value, index, self) {
	    return self.indexOf(value) === index;
	}
/*
	First, we get all checkpoints into an array and also set the current checkpoint
	as well as the previous checkpoint to 0. The checkpoints which are higher than the
	current checkpoint will be hidden, the chckpoints which are lower will be in the
	completed section.
*/
	var checkp = $(".checkpoint");
	var curCheckp = 0;
	var preCheckp = 0;

/*
	The core functions controlling what happens when a participant clicks the submit
	button on a checkpoint question. Currently few things happen, mainly the current
	section will be slided up and the next section will be slided down.
*/

	function checkSubmit() {
		//Current target (section that has been completed)
		var cTarget = checkp[curCheckp];
		curCheckp += 1
		//Next target (next checkpoint section)
		var nTarget = checkp[curCheckp];
		//Slide up all the children of the completed section
		$(cTarget).children().slideUp();
		//Get the data-cp attribute of the previous section and append it to the completed chapters
		$("#completed-chapters").append($("<h4 class='cp-header' style='cursor:pointer'>" + $(cTarget).attr("data-cp") + "</h4>").click(function(){
			//If it has already been clicked, slide it up
			if($(this).hasClass("is-clicked")){
				$(cTarget).children().slideUp();
				$(this).toggleClass("is-clicked");
			//Otherwise, slide it down
			} else {
				$(cTarget).children().slideDown();
				$(this).toggleClass("is-clicked");
			}
		}));
		//Slide down the next chapter
		$(nTarget).slideDown();
		//Reset the view to the highest point
		scrollTo(0, 0);
	}
/*
	The next block handles communication with socket.io. This is where the
	data gets sent to the node.js script and relayed to the chart page
*/
	//On form submit
	$('form').on("submit",function(f){
		//Clean the previous popup message (if any)
	//$("#answer-popup").empty();
		//Preventing default form behaviour like reloading
	//f.preventDefault();
	//f.stopPropagation();
		//This jQuery filter finds the active (checked) button(s) on the form
//		var activeBut = $(this).find("input").filter(":checked");
		//Setting up arays to hold the answers which are true as well as the
		//answer which were selectable
//		var TRUTH = [];
//		var TRUTHTEXTS = [];
//		var INPUTS = $(this).find("input");
//		var INPUTTEXTS = $(this).find(".i-text");
//		var COMMENT = $(this).find(".questcomment").text();
		//Iterate through all possible answers and add the truthful ones to the array
//		for(var i = 0; i < INPUTS.length; i++){
//				if($(INPUTS[i]).val()=="true"){
//					TRUTH.push(INPUTS[i]);
//					TRUTHTEXTS.push($(INPUTTEXTS[i]).text());
//				}
//			}
		//Iterate over each selected button in the form and send the data to the server
//		var i = 0;
//	    activeBut.each(function(){
		    //Setting up empty object
//		    var obj = {};
		    //To generate the question ID, remove all nasty special characters
//		    obj[$(this).attr("name").replace(/ /g, "\\ ").replace(/[\(\)\+\-\*\?\\\/]/g,'').replace(/[ \,\.]/g,'-')] = {
			    //The binary answer is the value of the input
//			    "binary":$(this).val(),
			    //The detailed answer is also stored in an attribute
//			    "detailed":$(this).attr("data-details"),
			    //We can get the cookie live
//			    "part":getPartCookie()
//		    }
		    //Send the data to the socket
//			socket.emit("quest answer",obj);
			//Since there is only one truthful answer per question we simply check
			//the detailed answer of the truth array vs. the one the participant chose
//			var ddat = $(activeBut[i]).attr("data-details");
			//If it matched send the correct popup
//			if(ddat === $(TRUTH[i]).attr("data-details")){
//				var curNam = $(TRUTH[i]).attr("name");
	//$("#answer-popup").append("<div class='corr-popup'> Thank you, close this popup to unfold the next section."+"</div>");
//			} else {
				//Otherwise fill the popup with the answer that would have been correct
//				$("#answer-popup").append("<div class='err-popup'> Q"+ (i+1) +": Oops! The correct answer was: "+ TRUTHTEXTS[i] + "</div>");
//			}
//			i++;
		//Also fade the background to gray while the popup is active
//	    });
//	    if(COMMENT != ""){
//			$("#answer-popup").append("<div class='questcomment'>" + COMMENT + "</div>");
//		}
    //$("#answer-popup-wrapper, #overlay").fadeIn();
		//Trigger animations
	    checkSubmit();
		//Also remove the button, to prevent multiple clicks, which would mess up the
		//section targeting...
	    $(this).find("button").remove();
	    return false;
	});
	//We can take all the questions and store them in an array
	var qs = [];
	//Clean up the names before pushing
	$("form").each(function(f){
		$(this).find("input").each(function(){
			qs.push($(this).attr("name").replace(/ /g, "\\ ").replace(/[\(\)\+\-\*\?\\\/]/g,''))
		})
	});
	//When initiating the participant, check if a cookie has been set, otherwise set it
	function initPart(callback){
		if(! checkPartCookie())
		{
                setPartCookie();
                callback();
    	}
    	}
    //Function to sent the participant hash to the server
    function emitParts(){
	    socket.emit("partconn",getPartCookie());
    }
    //Always send chart init message
    //FIXME:: Check if charts have been initiated, then send
    socket.emit("chart init",qs.filter(uniq));
	//Once the doc is loaded...
	$(document).ready(function(){
		//hide all checkpoints
		$(".checkpoint").hide();
		//un-hide the first checkpoint
		$(checkp[curCheckp]).show();
		//initiate the participant and send the hash
		initPart(emitParts);
		//send the hash again if the callback messed up
		emitParts();
		//Add the overlay div
		$("body").prepend("<div id='overlay'> </div>");
		//Fade out the divs that should not be visible
		$("#answer-popup-wrapper").click(function(){$(this).fadeOut(); $("#overlay").fadeOut()});
		$("body").prepend("<div id='partid'> </div>");
		$("#partid").text(getPartCookie());
	});
