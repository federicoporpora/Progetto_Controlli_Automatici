FILE README DA SCRIVERE PER RUNNARE GLI SCRIPT

Il file 'main.m' calcola i punti dal 3 al 4 e salva nella cartella ../report/figs i grafici necessari per compilare il pdf.
Per runnare lo script semplicemente aprire il file con matlab e premere run.
Per un risultato migliore nella visualizzazione del file pdf (se si vuole compilare) consigliamo di utilizzare il tema chiaro di matlab.

Il file 'simulink_lineare.slx' è la ricostruzione del sistema lineare con i disturbi. Per visualizzare il grafico finale, è necessario
prima runnare il file 'main.m' per caricare le variabili nell'ambiente, poi aprire il file 'simulink_lineare.slx' e runnarlo
con il solver ode45 ed una max step size di 1e-4 (dovrebbero essere già settate).
Il file ha il rumore disattivato per visualizzare in maniera chiara i grafici, ma è possibile inserire i dati richiesti dal testo per
verificare la correttezza del sistema anche in presenza di rumore.

Il file 'simulink_non_lineare.slx' è la ricostruzione del sistema non lineare con i disturbi.