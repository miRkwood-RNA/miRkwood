<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../Style/bioinfo.css' />
		<link type='text/css' rel='stylesheet' href='./style/script.css' />
		<link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
    
    <body>
        <div class="theme-border"></div>
        
        <a href="/">
            <div class="logo"></div>
        </a>
            
        <? include("./static/bioinfo_menu.txt") ?>
        
        <div class="bloc_droit">
            
            <? include("./static/header_menu.txt") ?>
            
            <div class="main">

                <br/>
                <h2>Retrieve result with an ID</h2>
                <br/>
  
  
                <div class="forms">
                    <form method="post" action="/cgi-bin/mirkwood/web_scripts/retrieve_job_with_ID.pl">
                        <p>The ID remains valid 15 days after sequence submission.</p>
                        <p>
                            <label for='jobId'><b>Enter the ID:</b></label>
                            <input type="text" name="jobId" id='jobId' size="20"/>
                            <input type="hidden" name="command" value="result"/>
                            <input type="submit" value="Go"/>
                        </p>
                    </form>
                </div> <!-- form -->


            </div> <!-- main -->
        </div><!-- bloc droit-->
        
        <? include("./static/footer.txt") ?>
       
    </body>
    
</html>
