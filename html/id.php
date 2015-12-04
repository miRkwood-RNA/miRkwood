<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>

        <link type='text/css' rel='stylesheet' href='../Style/css/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='./style/mirkwood.css' />
		<link type='text/css' rel='stylesheet' href='./style/script.css' />
		<link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type="text/javascript" src="../../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../../libs/jquery.history.js"></script>
        <script type="text/javascript" src="../../scripts/bioinfo.js"></script>

        <title>miRkwood - MicroRNA identification</title>
    </head>

    <body>
        <div class="frametitle">
            <h1 id="title">miRkwood</h1>
        </div>

        <div id="center_sup">
            <div id="link_home" style="display:inline-block"><a href="../index.php" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
            <div class="tabs" id="menu_central" style="display:inline-block">
                <? include("./static/header_menu.txt") ?>
            </div>
            <div id="arborescence"></div>
        </div>

        <div id="main">

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

        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>

    </body>

</html>
