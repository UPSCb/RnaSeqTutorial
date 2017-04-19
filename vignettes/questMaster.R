suppressMessages(require(shiny))
beginMainQuest <- function(){
  ret <- as.character(includeCSS("css/sock_style.css"))
  ret <- append(ret,"<div id='completed-chapters'></div>")
  #ret <- append(ret, "<div id='answer-popup-wrapper'><button>Close</button><div id='answer-popup'></div></div>")
  HTML(paste(ret, collapse=" "))
}
endMainQuest <- function(){
  includeScript("js/sock_fun.js")
}
startQuest <- function(title){
  HTML(paste("<div class=checkpoint data-cp='","+","    ",title,"'>"))
}
endQuest <- function(){
  HTML("</div>")
}
#quest <- function(l,lVec,comment=""){
quest <- function(n){
  ret <- character()
  #dVec <- c("A","B","C","D")
  ret <- append(ret,"<form>")
#  for(f in seq_along(l)){
#    tn <- names(l)[f]
    ret <- append(ret,paste("<hr/><br/><h3>",paste("Time for Socrative question number",n),"</h3><br/><hr/>",sep=""))
#    for( i in seq_along(l[[f]])){
#      val <- ifelse(lVec[f]==i,"true","false")
#      ret <- append(ret,paste('<input class="quest" type="radio" name="',tn,'" data-details="',dVec[i],'" value="',val,'">','<span class="i-text">',l[[f]][i],'</span>','<br/>',sep=""))
#    }
#  }
  ret <- append(ret,'<button class="s-btn">Click to proceed once you answered the question</button>')
  ret <- append(ret,'</form>')
                # paste(
                #   '<span class="questcomment" style="display: none">',
                #   comment,'</span></form>',sep = ""))
  HTML(ret)
}

