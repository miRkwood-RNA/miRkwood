/* These functions are taken from /bio1/www/html/scripts/bioinfo.js
 * miRkwood is not a dynamic website, so we cannot use the file
 * bioinfo.js (which is common to all softwares) because the anchors
 * on help pages wouldn't work anymore.
 */

var soft_theme = {
    /* theme DNA */
    "yass":"dna", 
    "magnolia":"dna", 
    "procars":"dna", 
    /* theme RNA */
    "carnac":"rna", 
    "gardenia":"rna", 
    "regliss":"rna", 
    "CGseq":"rna", 
    "mirkwood":"rna", 
    /* theme HTS */
    "sortmerna":"hts", 
    "crac":"hts", 
    "vidjil":"hts", 
    "storm":"hts", 
    /* theme PROTEINS */
    "path":"proteins", 
    "protea":"proteins", 
    "reblosum":"proteins", 
    /* theme TFM */
    "tfm-explorer":"tfm", 
    "tfm-scan":"tfm", 
    "tfm-pvalue":"tfm", 
    "tfm-cuda":"tfm", 
    /* theme NRP */
    "norine":"nrp", 
    "mynorine":"nrp", 
    "smiles2mono":"nrp", 
    "smilescolor":"nrp"
};


function is_a_soft_page (url) {
    var url_tab = url.split("/");
    return (url_tab[1] in soft_theme);
}

/*
 * Une fois que le DOM est chargé :
 * - vérifie si le site doit être dynamique ou non
 * - mets à jour le fil d'arianne si on est sur une page de logiciel
 * - active l'onglet
 */
 $(document).ready(function(){
    var url=window.location.pathname;
    // si c'est un logiciel, on affiche le fil d'arianne
    if (is_a_soft_page(url)) {
	display_arborescence(url);
    }
    set_appropriate_title(url);
    // active l'onglet
    setTimeout(function(){active_class(url);}, 200); // ??? pourquoi un timeout ???
});

/*
 * Affiche le fil d'arianne sur la page.
 *
 * Parametres
 *  url: URL complète d'une des pages d'un logiciel
 */
function display_arborescence(url){
    var tab_url=url.split("/");
    var name_soft=tab_url[1];
    if(soft_theme.hasOwnProperty(name_soft)){
		var theme = soft_theme[name_soft];
		$("#arborescence").html("<a href='/'>home</a> &gt; <a href='/theme_page/"+theme+".html'>"+theme.toUpperCase()+"</a> &gt; "+"<a href='/"+name_soft+"/"+name_soft.toLowerCase()+".php'>"+name_soft+"</a>");
    }
}

/*
 * Modifie le titre de la page dynamiquement;
 *
 * Parametres
 *  url: URL complète d'une des pages 
 */
function set_appropriate_title (url) {
    var tab_url = url.split("/");
    var name_soft = tab_url[1];
    var onglet = tab_url[2].replace(/.php/,"");
    if (is_a_soft_page (url)) {
        var title = "Bonsai bioinformatics - " + name_soft
        if (onglet != name_soft) {
            // on va chercher le texte de l'onglet affiché sur la page
            var elt_current=".tabs"; 
            var val = $(elt_current).find('a[href="'+url+'"]').text();
            title = title + " - " + val;
        }
        $(document).attr("title", title);
    } else if (name_soft == "theme_page") {
        onglet = onglet.replace(/.html/,"");
        var title = "Bonsai bioinformatics - " + onglet;
        $(document).attr("title", title);
    }	

}

/* 
 * Activer l'onglet sélectionne
 *
 * Parametres:
 *   url : l'URL nom du soft avec l'onglet a activer
 */
function active_class(url){
    //~ if (url.search("results.php") >= 0) {
        //~ var tab_url = url.split("/");
        //~ url = "/" + tab_url[1] + "/id.php";
    //~ }
    //~ alert("ACTIVE CLASS: " + url);
    // tabs est la classe de la div qui contient le menu central
    var elt_current=".tabs"; 
    // ajoute la classe active a l'element parent de l'element d'interet
    $(elt_current).find('a[href="'+url+'"]').parent().addClass("active");
    // supprime la classe active a tous les elements freres
    $(elt_current).find('a[href="'+url+'"]').parent().siblings().removeClass("active");
}
