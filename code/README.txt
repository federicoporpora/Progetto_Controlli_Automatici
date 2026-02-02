ISTRUZIONI PER RUNNARE GLI SCRIPT

FILE 'main.m' Questo è il file principale. Fa tutti i calcoli dei punti 3 e 4 e salva in automatico i grafici nella cartella ../report/figs per il pdf.
Come runnarlo: Apri con Matlab e premi Run.
Nota: Fatelo girare PRIMA di aprire i file Simulink, altrimenti mancano le variabili nel workspace e da errore.
Consiglio: Usate il tema chiaro di Matlab così i grafici escono giusti per il report.

FILE 'simulink_lineare.slx' E' lo schema del sistema lineare.
Prima di aprirlo: Assicuratevi di aver runnato 'main.m'.
Settings: Il file è già impostato con solver ode45, stop time a 0.5s e max step size a 1e-4 (serve per vedere bene la dinamica veloce).
Nota: Il rumore è disattivato di default per avere i grafici puliti, ma se serve si può collegare.

FILE 'simulink_non_lineare.slx' E' lo schema del sistema non lineare completo.
Prima di aprirlo: Sempre runnare 'main.m'.
Settings: Qui le impostazioni cambiano perché c'è il rumore ad alta frequenza. Lo stop time è a 2s e il max step size è 1e-6 (molto basso, altrimenti la simulazione sballa). Ci mette un attimo in più a finire rispetto all'altro.