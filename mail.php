<?php

$headers = "From: nadiroun@gmail.com <nadiroun@gmail.com>\n";
			$headers .= "MIME-Version: 1.0\n";
			$headers .= "Return-Path: <nadiroun@gmail.com>\n";
			$headers .= "Content-Type: text/html; charset=iso-8859-1\n";
			$headers .= "X-Sender: <monprojet.com>\n";
			$headers .= "X-Mailer: PHP\n";
			$headers .= "X-auth-smtp-user: contact@beloola.com\n";
			$to = "nadir@v-cult.com";
			$subject = "Demande inscription Beloola.com!";
			$body = "Bonjour, merci";
	

	
			$bodyretour = "<body><div style=\"text-align:justify;font-family: Verdana,Arial,Helvetica,sans-serif;\">
			Bonjour ,<br>Votre inscription &agrave; la version B&ecirc;ta-test de Beloola a bien &eacute;t&eacute; prise en compte avec le mot de passe. Votre compte est actuellement en cours de validation, vous receverez tr&agrave;s prochainement un mail vous informant de la validation de votre compte. Cette phase de B&ecirc;ta-test ferm&eacute;e permettra &agrave; nos &eacute;quipes, sur vos suggestions, de d&eacute;velopper au mieux les services propos&eacute;s au sein de Beloola. <br>Pour toute information ou question, contactez nous &agrave; l'adresse suivante : contact@beloola.com<br><br>
			<div align=\"center\"><img src=\"http://www.v-cult.com/upload/49.jpg\" alt=\"http://www.v-cult.com/upload/49.jpg\"></div>
			<br>...Et &agrave; bient&ocirc;t sur Beloola !</div>
			<br>
			<a href=\"http://www.beloola.com\" target=\"_blank\">www.beloola.com</a><br>
			<a href=\"http://www.beloola.com/blog\" target=\"_blank\">www.beloola.com/blog</a>

			<span style=\"font-size: 13.3px; font-family: Verdana,Arial,Helvetica,sans-serif;\"><div style=\"margin: 0pt 0pt 8px;\">
			</div><span style=\"color: gray;\">&nbsp;</span><a target=\"_blank\" style=\"text-decoration: underline;\" href=\"http://www.facebook.com/beloola\"><img border=\"0\" width=\"16\" height=\"16\" src=\"http://images.wisestamp.com/facebook.png\" style=\"padding: 0px 0px 5px; vertical-align: middle;\" alt=\"Facebook\"></a> <a target=\"_blank\" style=\"text-decoration: underline;\" href=\"http://www.twitter.com/beloola\"><img border=\"0\" width=\"16\" height=\"16\" src=\"http://images.wisestamp.com/twitter.png\" style=\"padding: 0px 0px 5px; vertical-align: middle;\" alt=\"Twitter\"></a> <a target=\"_blank\" style=\"text-decoration: underline;\" href=\"http://www.beloola.com/blog/rss.php\"><img border=\"0\" width=\"16\" height=\"16\" src=\"http://images.wisestamp.com/blogRSS.png\" style=\"padding: 0px 0px 5px; vertical-align: middle;\" alt=\"Blog RSS\"></a><br>
			</span></body>";

		var_dump(mail("nadiroun@gmail.com", "Inscription Beloola", $bodyretour,$headers));
			
	?>
