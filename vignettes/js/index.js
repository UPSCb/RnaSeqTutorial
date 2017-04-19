var app = require('express')();
var http = require('http').Server(app);
var io = require('socket.io')(http);

app.get('/', function(req, res){
  res.sendFile(__dirname + '../index.html');
});

http.listen(3000, function(){
  console.log('listening on *:3000');
});

/*
Define a callback function to filter out redundant values of an array
in order to only keep unique participant hashes.
*/
function uniq(value, index, self) { 
	    return self.indexOf(value) === index;
	}

/*
Backup variable to store transmitted data on the server and echo back 
on initialisation request.	
*/
var syncdata = {"partconn" : [], "chdat" : [], "chinit" : []};

/*
Here, actions that happen upon socket connections are defined. Each
action gets an ID like "quest answer" to indicate what the connection
is meant for.	
*/
io.on('connection', function(socket){	
	/*
	This function echoes back all synchronized data on all pages upon
	an initialisation request from the chart page (once). Currently
	badly implemented cause no check to remove redundancy in the 
	participant answer data (will still only count as one in the
	final page, but waste of resources).	
	*/
	socket.on('init request', function(){
		/*
		Echo once for participant connection data and once for the
		question answer data to plot the charts.
		*/
		for(var i = 0; i < syncdata['partconn'].length; i++){
			io.emit('partconn', syncdata['partconn'][i]);
		}
		for(var i = 0; i < syncdata['chinit'].length; i++){
			io.emit('chart init', syncdata['chinit'][i]);
		}
		for(var i = 0; i < syncdata['chdat'].length; i++){
			io.emit('quest answer', syncdata['chdat'][i]);
		}
	});
	/*
	On question answer by the participant, log the data on the server
	and then simply echo the answer data to the chart page.
	Current format of `msg`: {"Question" : {"binary":,"detailed":,"part"}},
	where `binary` stores true or false, `detailed` stores A, B, C or D and
	`part` stores the participant's hash.
	*/
	socket.on('quest answer', function(msg){
		console.log(msg);
		var k = Object.keys(msg);
		syncdata['chdat'].push(msg);
	    io.emit('quest answer', msg);
	});
	/*
	On participant connection, simply echo the hash back to the chart server.	
	*/
	socket.on('partconn', function(msg){
		syncdata['partconn'].push(msg);
		//filter redundant hashes
		syncdata['partconn'] = syncdata['partconn'].filter(uniq);
    	io.emit('partconn', msg);
    });
    /*
	Auto initialisation of all charts by communication with the client.
	*/
	socket.on('chart init',function(msg){
		syncdata['chinit'].push(msg);
		syncdata['chinit'] = syncdata['chinit'].filter(uniq);
		io.emit('chart init',msg);
	});
    
});