# Projet-DSP-C6713

Ce projet permet de filtrer un signal. On estime sur des trames composés uniquement de bruit la puissance spectrale moyenne du bruit. On décompose le signal d'entrée en trames de 256 échantillons dans lesquels on soustrait la puissance spectrale moyenne du bruit(pondéré d'un facteur) à la puissance spectrale des trames. Pour supprimer le bruit spectral qui pourrait apparaitre, on réinsère du bruit.

Ce filtrage est codé en matlab et en C.
Le code C fonctionne (entre autre) sur le DSP de Texas Instrument : C6713.

Le code utilise des fonctions provenant d'une librairie de TI.

Le code n'est pas optimisé ni très propre mais fonctionne.
