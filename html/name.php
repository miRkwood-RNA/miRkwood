<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='./style/mirkwood.css' />
       <!-- <link type='text/css' rel='stylesheet' href='./style/help.css' /> -->
        <link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <link type='text/css' rel='stylesheet' href='../Style/css/bioinfo.css' />
        <script type="text/javascript" src="../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="../libs/jquery.history.js"></script>
        <script type="text/javascript" src="./js/header.js"></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
    <body>

        <div class="frametitle">
            <h1 id="title">miRkwood</h1>
        </div>

        <div class="theme-border" style="display:none"></div>
        <div class="tabs" id="menu_central" style="display:inline-block"><?php include("./static/header_menu.txt")?></div>
        <div id="arborescence"></div>

        <div id="main">

            <h2>Where does the name come from ?</h2>

            <br /><br />

            <img style='width:600px; display: block; margin: 0 auto;'src='./style/forest.jpg' alt='mirkwood forest' />

            <br /><br />

            <p>The term Mirkwood is taken from the germanic mythology (Myrkvior, meaning "mirky wood, dark wood")
            and is used for two distinct fictional forests on the continent of Middle-earth in Tolkien's legendarium.
            </p>

            <p>One of these occurred in the First Age of Middle-earth, when the highlands of Dorthonion north
            of Beleriand were known as Mirkwood after falling under Morgoth's control.
            The other Mirkwood, and the more famous of the two, was the large forest in Rhovanion, east of the Anduin.
            This had acquired the name Mirkwood during the Third Age, after it fell under the influence of the Necromancer.
            This Mirkwood features significantly in <i>The Hobbit</i> and in the film <i>The Hobbit: The Desolation of Smaug</i>.
            </p>

        </div>

        <?php include('static/footer.php') ?>

    </body>

</html>
