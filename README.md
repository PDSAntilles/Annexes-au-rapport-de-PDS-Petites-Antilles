# Annexes - Rapport de PDS: Détection de séismes répétés dans la subduction des Petites Antilles par template-matching (2025)#

> Ce dépôt contient 6 fichiers, dont le rapport final de PDS en .pdf, 3 fichiers de code annexes reliés à l'étude présentée dans le rapport, et 2 grandes images plus simplement visualisables sur cette page Git.
> Les trois fichiers sont des versions modifiées des codes de C. Satriano sur Requake (https://github.com/SeismicSource/requake).
> Bien que la majorité des fonctionnalités de Requake aient été utilisées telles quelles, certains fichiers ont été modifiés pour adapter 
  les calculs de template_matching aux longues périodes étudiées dans le cadre de ce PDS. Les codes utilisés sans modifications sont accessibles via le précédent lien GitHub.
> Les images sont accessibles directement au format PNG.

# Contenu du dépôt:

  /docs

     1/ scan_templates.py 
  
     2/ waveforms.py
  
     3/ requake.conf

  /Images
  
     1/ famille-.png

     2/ famille+.png

  /Rapport_de_stage_IPGP_LE_CAM_YVARD_2025.pdf

     
# Description du dépôt:

 Fichier docs/1 : Ce fichier .py permet via la fonction scan_templates() de lancer un calcul de recherche par templates sur une longue plage de temps définie dans requake.conf.
     Il utilise un template au format .sac, associé à une famille supposée. La chaîne de calcul associée est détaillée dans la section 2.3 du rapport PDF.

 Fichier docs/2 : Ce second fichier .py gère les requêtes aux serveurs afin de renvoyer des données sismiques continues, exploitables par scan_templates.py . 
     Le fichier gère aussi l'alignement des traces et les calculs mathématiques précédent l'affichage des traces (notamment les templates moyens) 
     via la commande "requake plot_families".

 Fichier docs/3 : Ce fichier .conf gère toutes les variables de type float ou string permettant de paramétrer les calculs de Requake.

 Images Image/famille-.png : Ensemble des 20 évènements de la famille dite "famille -" (Voir section 3.3 du rapport)

 Images Image/famille+.png : Ensemble des 30+3 évènements de la famille dite "famille +" (Voir section 3.3 du rapport), avant réduction à 30 évènements

      ---> Les deux images sont directement associées à la figure 40 du rapport, où le déplacement local de la source de     repeaters de la famille considérée est encadré par deux valeurs : une borne inférieure (famille-), et une borne supérieure (famille+)

  Fichier /Rapport_de_stage_IPGP_LE_CAM_YVARD_2025.pdf : Rapport final de PDS (50 p.)
 



