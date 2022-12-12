\NeedsTeXFormat{LaTeX2e}

\RequirePackage[polish]{babel}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}

\def\@baseclass{report}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\@baseclass}}
\PassOptionsToClass{a4paper, twoside, 11pt}{\@baseclass} %openright
\ProcessOptions

\LoadClass{\@baseclass}

\textwidth
\paperwidth
\advance\textwidth -55mm
\oddsidemargin-1in
\advance\oddsidemargin 30mm
\evensidemargin-1in
\advance\evensidemargin 25mm
\topmargin -1in
\advance\topmargin 2cm
\setlength\textheight{48\baselineskip}
\addtolength\textheight{\topskip}
\marginparwidth15mm


\renewcommand
\maketitle
{%
	\begin{titlepage}%
		\begin{center}%
			{
				\LARGE\textbf{Uniwersytet Warszawski}\\
				Wydzia\l{} Fizyki\par
			}
			\vspace{1cm plus .5fill}
			{
				\centering
				\begin{minipage}{0.45\textwidth}
					\centering{\Large\bfseries Oleg Kmechak}\\
					[2mm]Nr albumu: 383501
				\end{minipage}
			}

			\vspace{8mm plus .5fill}
			{
				\Huge
				\textbf{Modeling of the growth of river networks using the Riversim package}
				\par
			}

			\vspace{8mm plus .1fill}
			{
				\large\bf Praca magisterska\\[3pt]
				na kierunku FIZYKA \\
			}
			\vspace{2cm plus 1.5fill}
			\begin{flushright}\large
				\begin{tabular}{l}
					Praca wykonana pod kierunkiem\\
					\bfseries prof. dr hab. Piotr Szymczak\\
					Instytut Fizyki Teoretycznej\\
				\end{tabular}
			\end{flushright}
			\vspace{1cm plus 1fill}
			{
				\large Warszawa, Wrzesień 2022\par
			}
		\end{center}
	\end{titlepage}
}


\renewenvironment{abstract}
{
	\begin{center} \bfseries \large \abstractname \end{center}
}
{
	\par
	\vspace*{26pt}%
	\begin{center} \bfseries\large Słowa kluczowe \end{center}
  	FEM, Fizyka, Wody podziemne, Rownanie Laplace'a, Erozja, Metody numeryczne\par
  	
	\vspace*{26pt}%
	\begin{center} \bfseries\large Tytuł pracy w języku angielskim \end{center}
	Modelowanie wzrostu sieci rzecznych za pomocą pakietu Riversim.\par
}