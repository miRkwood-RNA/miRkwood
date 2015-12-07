<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../Style/css/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='./style/mirkwood.css' />
        <link type='text/css' rel='stylesheet' href='./style/help.css' />
        <link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type="text/javascript" src="../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../libs/jquery.history.js"></script>
        <script type="text/javascript" src="./js/header.js"></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
    <body>

        <div class="frametitle">
            <h1 id="title">miRkwood</h1>
        </div>

        <div id="center_sup">
            <div class="theme-border" style="display:none"></div>
            <div id="link_home" style="display:inline-block"><a href="/theme_page/rna.html" class="text_onglet"><img src="/Style/icon/home_w.png" alt="home_general"/></a></div>
            <div class="tabs" id="menu_central" style="display:inline-block"><?php include("./static/header_menu.txt")?></div>
            <div id="arborescence"></div>
        </div>

        <div id="main">

            <p class='title'>Where does the name come from ?</p>

            <img style='width:600px; display: block; margin: 0 auto;'src='./style/forest.jpg' alt='mirkwood forest' />

            <br />

            <p>
            The term Mirkwood is taken from the germanic mythology (Myrkvior, meaning "mirky wood, dark wood") and is used  for two distinct fictional forests on the continent of Middle-earth in Tolkien's legendarium.
            </p>

            <p>
            One of these occurred in the First Age of Middle-earth, when the highlands of Dorthonion north of Beleriand were known as Mirkwood after falling under Morgoth's control. The other Mirkwood, and the more famous of the two, was the large forest in Rhovanion, east of the Anduin. This had acquired the name Mirkwood during the Third Age, after it fell under the influence of the Necromancer. This Mirkwood features significantly in The Hobbit and in the film The Hobbit: The Desolation of Smaug.
            </p>

        </div>

        <?php require("/bio1/www/html/lib.inc")?>
        <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>

    </body>

</html>
