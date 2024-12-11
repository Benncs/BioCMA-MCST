
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{docclass}[2024/01/09 Class for package documentation]

% Load standard article class with options
\LoadClass[11pt, a4paper]{article}
\newcommand{\myversion}{0.3.x}


% Load necessary packages
\usepackage[margin=0.6in]{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{amsmath}
\usepackage{mathtools}
\usepackage{float}
\usepackage{stmaryrd}
\usepackage[normalem]{ulem}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage[backend=bibtex,style=authoryear,maxnames=2,natbib=true]{biblatex}
\usepackage{gensymb} %Degree C 
\usepackage{fancyhdr}
\usepackage{amsfonts}
\usepackage{bm}
\numberwithin{equation}{section}

% Define figure command
\newcommand{\mkfig}[3]{%
    \begin{figure}[h]
        \centering
        \includegraphics[width=0.8\textwidth, keepaspectratio]{#1}
        \caption{\string#2}
        \label{fig:#3}
    \end{figure}
}

% Define citation figure command
\newcommand{\lcitefig}[3]{%
   \hyperref[#1]{\string#2 \ref{#3}}
}

\definecolor{backcolour}{rgb}{0.95,0.95,0.92}
\definecolor{codegreen}{rgb}{0,0.6,0}
\definecolor{codegray}{rgb}{0.5,0.5,0.5}
\definecolor{codepurple}{rgb}{0.58,0,0.82}
\definecolor{backcolour}{rgb}{0.95,0.95,0.92}


\lstdefinestyle{code_style}{
    backgroundcolor=\color{backcolour},
    commentstyle=\color{codegreen},
    keywordstyle=\color{codepurple},
    numberstyle=\tiny\color{codegray},
    stringstyle=\color{codepurple},
    basicstyle=\ttfamily\footnotesize, % Smaller font size
    breakatwhitespace=false,
    breaklines=true,
    captionpos=b,
    keepspaces=true,
    showspaces=false,
    showstringspaces=false,
    showtabs=false,
    tabsize=2,
    framexleftmargin=0pt, % Full width
    frame=shadowbox,
    rulesepcolor=\color{gray},
    xleftmargin=0pt,
    frame=single,
    framextopmargin=1mm, % Smaller top margin
    framexbottommargin=1mm, % Smaller bottom margin
    framesep=1pt,
    linewidth=\textwidth
}

\newcommand{\addListofFiguresTablesBibliography}{
  \clearpage 
  \addcontentsline{toc}{section}{\listfigurename}
  \listoffigures
  \addcontentsline{toc}{section}{\listtablename}
  \listoftables
  \addcontentsline{toc}{section}{References}
  \printbibliography
}

% Set code listing style
\lstset{style=code_style}

% Load biblatex with references file
\addbibresource{../common_ref}

\addbibresource{references}



\pagestyle{fancy}
\fancyhf{} % Clear default headers and footers
\rfoot{Version \myversion} % Right-aligned footer

