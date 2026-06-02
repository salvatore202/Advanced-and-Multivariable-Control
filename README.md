# Advanced-And-Multivariable-Control

[![MATLAB/Simulink](https://img.shields.io/badge/MATLAB-Simulink-orange.svg)](https://it.mathworks.com/products/simulink.html)
[![Platform-Web](https://img.shields.io/badge/Platform-Web%20%2F%20HTML5-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

Raccolta delle prove d'esame svolte e del materiale fornito dal docente del corso di **Complementi di Controlli- Controllo Avanzato e Applicazioni (6 CFU)**. Ogni directory del tipo Esame_*codice* include:
1. file script Main_*codice*  `.m` con il codice matlab
2. file live script Main_*codice*  `.mlx`
3. file Main_*codice*  `.ipynb` con il live script consultabile direttamente su questa pagna web
4. file Simulink Verifica_*TYPE*_*codice*  `.slx`, con le gli schemi simulink di verifica di funzionamento del controllore
5. file Traccia_*codice* `.pdf`
6. **Opzionalmente:** file immagini `.png` raffiguranti immagini per la compressione SVD o schemi di controllo della traccia



## Struttura della Repository

La repository è organizzata in 4 directory principali. Indicazioni su struttura e contenuto delle directory di seguito
```text
├── .github/                     # Configurazioni per GitHub Pages
│
├── Esame/                       # Directory usata per creare il file di consegna 
│
├── Esercizi_matlab/             # Esercitazioni e Prove d'esame svolte
│   ├── Esame_11032026/
│   ├── Esame_11052026/
│   ├── Esame_20022026/
│   ├── Esame_20042026/
│   ├── Esame_27042026/
│   └── Esercitazioni/
│
├── Materiale_Teams/            # Materiale condiviso dal prof. durante il corso
│   ├── 1. STABILITA'/
│   ├── 2. STABILITA' ROBUSTA/
│   ├── 3. OSSERVABILITA' E RAGGIUNGIBILITA'/
│   ├── 4. STABILIZZAZIONE E ASSEGNAMENTO AUTOVALORI/
│   ├── 5. SVD/
│   ├── 6. STABILITA' E STABILITà ROBUSTA NEL DOMINIO DI LAPLACE/
│   ├── 7. CONTROLLO Hinf/
│   ├── 8. CONTROLLO_OTTIMO/
│   └── 9. ESERCITAZIONE/
│
├── docs/                       # Report Web Schemi simulink (tutorial di seguito)
│
└── README.md                   # Questo file
```



## 🌐 Come Consultare i Report Interattivi (GitHub Pages)

Per ogni prova d'esame è stata generata una **Simulink Web View** per ogni file di Verifica Simulink in modo da consultare lo schema simulink direttamente da questa pagina web. Grazie a GitHub Pages, puoi esplorare gli schemi a blocchi Simulink e i risultati direttamente dal tuo browser seguendo il link di seguito, **senza bisogno di avere MATLAB installato**. 

https://salvatore202.github.io/Advanced-and-Multivariable-Control/



### Video dimostrativo

<img src="docs/tutorial_1.gif" alt="Come navigare nei modelli Simulink via Web" width="100%" />

## 📌 Guida Passo-Passo per l'utilizzzo della repository

#### Opzione 1: Visualizzazione Online (Consigliato)
1. Spostarsi in Esame_*codice*/
2. Aprire il file Main_*codice* `.m` per consultare il codice "row"
3. Aprire il file Main_*codice* `.ipynb` per consultare il live script della prova d'esame direttamente dal browser
4. Seguire il tutorial precedente per visualizzare e analizzare gli schemi simulink. **Attenzione:** le Web pages degli schemi simulink permettono solo di visualizzare la struttura dello schema simulink e dei blocchi. Non sarà possibile simulare (seguire Opzione 2 per effettuare simulazioni)

#### Opzione 2: Clone Locale e Riproduzione Offline (per Edit e Simulazioni)
Se desideri scaricare il materiale per eseguire gli script `.m`, i live scripts `.mlx` o modificare i modelli `.slx` in locale su MATLAB:

1. Apri il terminale e clona la repository:

   ```bash
   git clone https://github.com/salvatore202/Advanced-And-Multivariable-Control.git
   cd Advanced-And-Multivariable-Control
2. Avvia Matlab
3. **Installa dipendenze**:  cvx -> https://cvxr.com/cvx/download/
4. Nel caso in cui sia necessario usare cvx per eseguire lo script, scrivere nella Command Window di Matlab:

    ```bash
    cd ~/*percorso-cvx*/
    cvx_setup
    ```
5. Premere RUN nella sezione Live Editor di Matlab
6. **SI CONSIGLIA LA VISUALIZZAZIONE DEGLI OUTPUT INLINE**