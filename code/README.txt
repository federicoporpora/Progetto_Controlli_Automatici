ISTRUZIONI PER RUNNARE GLI SCRIPT

FILE 'main.m' File principale. Esegue i punti 3 e 4 e salva in automatico i grafici nella cartella ../report/figs per il pdf.
Per runnarlo aprire con Matlab e premere Run.
Da runnare prima dei vari file Simulink, altrimenti mancano le variabili nel workspace e da errore.
Consigliamo l'utilizzo del tema chiaro di Matlab così i grafici escono giusti per il report.

FILE 'simulink_lineare.slx' E' lo schema del sistema lineare.
Prima di aprirlo assicurarsi di aver runnato 'main.m'.
Il file dovrebbe essere già impostato con solver ode45, stop time a 0.5s e max step size a 1e-4 (serve per vedere bene la dinamica veloce).

FILE 'simulink_non_lineare.slx' E' lo schema del sistema non lineare completo.
Prima di aprirlo assicurarsi di aver runnato 'main.m'.
Qui le impostazioni cambiano perché c'è il rumore ad alta frequenza. Lo stop time è a 2s e il max step size è 1e-6 (molto basso per una maggiore precisione).
