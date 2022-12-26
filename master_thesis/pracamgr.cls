\NeedsTeXFormat{LaTeX2e}

\RequirePackage[english]{babel}
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
				\LARGE\textbf{University of Warsaw}\\
				Faculty of Physics\par
			}
			\vspace{1cm plus .5fill}
			{
				\centering
				\begin{minipage}{0.45\textwidth}
					\centering{\Large\bfseries Oleg Kmechak}\\
					[2mm]Record book number: 383501
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
				\large\bf Master's thesis\\[3pt]
				in the field of Physics \\
			}
			\vspace{2cm plus 1.5fill}
			\begin{flushright}\large
				\begin{tabular}{l}
					The thesis was written under the supervision of\\
					\bfseries prof. dr hab. Piotr Szymczak\\
					Chair of Complex Systems Modeling\\
					Institute of Theoretical Physics
				\end{tabular}
			\end{flushright}
			\vspace{1cm plus 1fill}
			{
				\large Warsaw, December 2022\par
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
	\begin{center} \bfseries\large Keyword \end{center}
  	FEM, physics, groundwater dynamics, erosion, numerical methods, network geometry\par
  	
	\vspace*{26pt}%
	\begin{center} \bfseries\large Title of the thesis in Polish language\end{center}
	Modelowanie wzrostu sieci rzecznych za pomocą pakietu Riversim.\par
}