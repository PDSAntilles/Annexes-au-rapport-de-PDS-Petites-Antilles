# Annexes - Rapport de PDS: Détection de séismes répétés dans la subduction des Petites Antilles par template-matching (2025)#

> Ce dépôt contient 3 fichiers annexes reliés à l'étude présentée dans le rapport.
> Ces trois fichiers sont des versions modifiées des codes de C. Satriano sur Requake (https://github.com/SeismicSource/requake).
> Bien que la majorité des fonctionnalités de Requake aient été utilisées telles quelles, certains fichiers ont été modifiés pour adapter 
  les calculs de template_matching aux longues périodes étudiées dans le cadre de ce PDS. Les codes utilisés sans modifications sont accessibles via le précédent lien GitHub.

# Contenu du dépôt:

  1/ scan_templates.py 
  
  2/ waveforms.py
  
  3/ requake.conf

# Description du dépôt:

 Fichier 1/ : Ce fichier .py permet via la fonction scan_templates() de lancer un calcul de recherche par templates sur une longue plage de temps définie dans requake.conf.
     Il utilise un template au format .sac, associé à une famille supposée. La chaîne de calcul associée est détaillée dans la section 2.3 du rapport PDF.

 Fichier 2/ : Ce second fichier .py gère les requêtes aux serveurs afin de renvoyer des données sismiques continues, exploitables par scan_templates.py . 
     Le fichier gère aussi l'alignement des traces et les calculs mathématiques précédent l'affichage des traces (notamment les templates moyens) 
     via la commande "requake plot_families".

 Fichier 3/ : Ce fichier .conf gère toutes les variables de type float ou string permettant de paramétrer les calculs de Requake.



