\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mdpi}[2017/11/21 MDPI paper class]
 
%% PRELIMINARY DECLARATIONS
\LoadClass[10pt]{article}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{calc}
\RequirePackage{indentfirst}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx,epstopdf}
\RequirePackage{lastpage}
\RequirePackage{ifthen}
\RequirePackage{lineno}
\RequirePackage{float}
\RequirePackage{amsmath}
\RequirePackage{setspace}
\RequirePackage{enumitem}
\RequirePackage{mathpazo}
\RequirePackage{booktabs} % For \toprule etc. in tables
\RequirePackage[largestsep]{titlesec}
\RequirePackage{etoolbox} % For \AtBeginDocument etc.
\RequirePackage{tabto} % To use tab for alignment on first page
\RequirePackage{color} % To provide color for soul
\RequirePackage{soul} % To highlight text
\newcommand{\highlight}[1]{\colorbox{yellow}{#1}}
\RequirePackage{multirow}
\RequirePackage{microtype} % For command \textls[]{}
\RequirePackage{tikz} % For Orcid icon

%% To choose the whether there is one or more authors
\DeclareOption{oneauthor}{\gdef\@authornum{author}}
\DeclareOption{moreauthors}{\gdef\@authornum{authors}}

%% Defaults
\ExecuteOptions{10pt,a4paper,article,oneauthor}

%% Process options
\ProcessOptions\relax

%% MORE DECLARATIONS
%%%% Maths environments
\RequirePackage{amsthm}
\newtheoremstyle{mdpi}% name
{12pt}% space above
{12pt}% space below
{\itshape}% body font
{}% indent amount 1
{\bfseries}% theorem head font
{.}% punctuation after theorem head
{.5em}% space after theorem head
{}% theorem head spec (can be left empty, meaning `normal')

\renewcommand{\qed}{\unskip\nobreak\quad\qedsymbol} %% This places the symbol right after the text instead of placing it at the end on the line.

\renewenvironment{proof}[1][\proofname]{\par %% \proofname allows to have "Proof of my theorem"
  \pushQED{\qed}%
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        \bfseries %% "Proof" is bold
    #1\@addpunct{.}]\ignorespaces\footnotesize %% Period instead of colon 
}{%
  \normalsize\popQED\endtrivlist\@endpefalse
}

 \theoremstyle{mdpi}
 \newcounter{theorem}
 \setcounter{theorem}{0}
 \newtheorem{Theorem}[theorem]{Theorem}
 
 \newcounter{lemma}
 \setcounter{lemma}{0}
 \newtheorem{Lemma}[lemma]{Lemma}
 
 \newcounter{corollary}
 \setcounter{corollary}{0}
 \newtheorem{Corollary}[corollary]{Corollary}
 
 \newcounter{proposition}
 \setcounter{proposition}{0}
 \newtheorem{Proposition}[proposition]{Proposition}
 
 \newcounter{problem}
 \setcounter{problem}{0}
 \newtheorem{Problem}[problem]{Problem}
 
 \newcounter{example}
 \setcounter{example}{0}
 \newtheorem{Example}[example]{Example}
 
 \newcounter{remark}
 \setcounter{remark}{0}
 \newtheorem{Remark}[remark]{Remark}
 
 \newcounter{definition}
 \setcounter{definition}{0}
 \newtheorem{Definition}[definition]{Definition}
 
 \newcounter{hypothesis}
 \setcounter{hypothesis}{0}
 \newtheorem{Hypothesis}[hypothesis]{Hypothesis}

 \newcounter{notation}
 \setcounter{notation}{0}
 \newtheorem{Notation}[notation]{Notation}
 
%%%% Hyphenation
\RequirePackage[none]{hyphenat}
\sloppy

%%%% References
\RequirePackage[sort&compress]{natbib}\bibliographystyle{chicago2}
\bibpunct{(}{)}{;}{x}{}{}
\renewcommand\NAT@set@cites{%
  \ifNAT@numbers
    \ifNAT@super \let\@cite\NAT@citesuper
       \def\NAT@mbox##1{\unskip\nobreak\textsuperscript{##1}}%
       \let\citeyearpar=\citeyear
       \let\NAT@space\relax
       \def\NAT@super@kern{\kern\p@}%
    \else
       \let\NAT@mbox=\mbox
       \let\@cite\NAT@citenum
       \let\NAT@space\relax
       \let\NAT@super@kern\relax
    \fi
    \let\@citex\NAT@citexnum
    \let\@biblabel\NAT@biblabelnum
    \let\@bibsetup\NAT@bibsetnum
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@num\NAT@close}%
    \def\natexlab##1{}%
    \def\NAT@penalty{\penalty\@m}%
  \else
    \let\@cite\NAT@cite
    \let\@citex\NAT@citex
    \let\@biblabel\NAT@biblabel
    \let\@bibsetup\NAT@bibsetup
    \let\NAT@space\NAT@spacechar
    \let\NAT@penalty\@empty
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@date\NAT@close}%
    \def\natexlab##1{##1}%
  \fi}

%%%%% Hyperlinks
%% Define color for citations
\definecolor{bluecite}{HTML}{0875b7}

\RequirePackage[unicode=true,
	bookmarksopen={true},
	pdffitwindow=true, 
	colorlinks=true, 
	linkcolor=bluecite, 
	citecolor=bluecite, 
	urlcolor=bluecite, 
	hyperfootnotes=false, 
	pdfstartview={FitH},
	pdfpagemode= UseNone
	]{hyperref}


%% To have the possibility to change the urlcolor
\newcommand{\changeurlcolor}[1]{\hypersetup{urlcolor=#1}} 

%% Metadata
\newcommand{\org@maketitle}{}% LATEX-Check
\let\org@maketitle\maketitle
\def\maketitle{%
	\hypersetup{
		pdftitle={\@Title},
		pdfsubject={\@abstract},
		pdfauthor={\@AuthorNames}
		}%
	\org@maketitle
}

%%%% Footnotes
\RequirePackage[hang]{footmisc}
\setlength{\skip\footins}{1.2cm}
\setlength{\footnotemargin}{5mm}
\def\footnoterule{\kern-14\p@
\hrule \@width 2in \kern 11.6\p@}

%%%% URL
\RequirePackage{url}
\urlstyle{same}
\g@addto@macro{\UrlBreaks}{\UrlOrds} 

%%%% Widows & orphans
\clubpenalty=10000
\widowpenalty=10000
\displaywidowpenalty=10000

%%%% Front matter
\newcommand{\firstargument}{}
\newcommand{\Title}[1]{\gdef\@Title{#1}}%
\newcommand{\Author}[1]{\gdef\@Author{#1}}%
\def\@AuthorNames{}
\newcommand{\AuthorNames}[1]{\gdef\@AuthorNames{#1}}%
\newcommand{\firstpage}[1]{\gdef\@firstpage{#1}}
\newcommand{\address}[2][]{\renewcommand{\firstargument}{#1}\gdef\@address{#2}}
\newcommand{\corresfirstargument}{}
\def\@corres{}
\newcommand{\corres}[2][]{\renewcommand{\corresfirstargument}{#1}\gdef\@corres{#2}} 
\def\@abstract{}
\renewcommand{\abstract}[1]{\noindent \textbf{Abstract. } \@abstract{#1}}




\def\@firstnote{}
\newcommand{\firstnote}[1]{\gdef\@firstnote{#1}}
\def\@secondnote{}
\newcommand{\secondnote}[1]{\gdef\@secondnote{#1}}%
\def\@thirdnote{}
\newcommand{\thirdnote}[1]{\gdef\@thirdnote{#1}}%
\def\@fourthnote{}
\newcommand{\fourthnote}[1]{\gdef\@fourthnote{#1}}%
\def\@fifthnote{}
\newcommand{\fifthnote}[1]{\gdef\@fifthnote{#1}}%
\def\@sixthnote{}
\newcommand{\sixthnote}[1]{\gdef\@sixthnote{#1}}%
\def\@seventhnote{}
\newcommand{\seventhnote}[1]{\gdef\@seventhnote{#1}}%
\def\@eighthnote{}
\newcommand{\eighthnote}[1]{\gdef\@eighthnote{#1}}%

\def\@simplesumm{}
\newcommand{\simplesumm}[1]{\gdef\@simplesumm{#1}}

\def\@externalbibliography{}
\newcommand{\externalbibliography}[1]{\gdef\@externalbibliography{#1}}

\def\@reftitle{}
\newcommand{\reftitle}[1]{\gdef\@reftitle{#1}}



%%%% Maketitle part 1: Logo, Arttype, Title, Author
\renewcommand{\@maketitle}{
	\begin{flushleft}
   	 	\vspace{-1pt}
  	  	\fontsize{18}{18}\selectfont
   	 	\centering\boldmath\bfseries{\@Title}
   	 	\par
   	 	\vspace{15pt}
   	 	 \fontsize{10}{10}\selectfont
   		\boldmath\bfseries{\@Author}
    	\par
    	\vspace{-4pt}
	\end{flushleft}%
	}

% Commands for hanging indent
\newcommand{\dist}{1.7em}
\newcommand{\hang}{\hangafter=1\hangindent=\dist\noindent}

%%%% Maketitle part 2
\newcommand{\maketitlen}{ 
	\begin{flushleft}
	\begin{spacing}{1.35}
    \vspace{15pt}
	\leftskip0.2cm
	\fontsize{9}{9}\selectfont
    \@address
    \hang\textbf{*} \tabto{\dist} \@corres
    \vspace{10pt}
    \par
    Version: {\@ \today}
	\end{spacing}
	\end{flushleft}
}


%%%% Abstract, keywords, journal data, PACS, MSC, JEL
%%%% Print maketitle and abstractkeywords
	\AfterEndPreamble{
	\maketitle
	\let\maketitle\relax
	\maketitle
	\maketitlen
	}%

\AtBeginDocument{
	\DeclareSymbolFont{AMSb}{U}{msb}{m}{n}
	\DeclareSymbolFontAlphabet{\mathbb}{AMSb}
	}

%%%% Font size in Tables
\AtEndPreamble{
	\def\@tablesize{}
	\newcommand{\tablesize}[1]{\gdef\@tablesize{#1}}
	\let\oldtabular\tabular
	\renewcommand{\tabular}{\ifthenelse{\equal{\@tablesize}{\@empty}}{\tiny}{\@tablesize}\oldtabular}
}

%%%% Section headings
\setcounter{secnumdepth}{4} %i.e., section numbering depth, which defaults to 3 in the article class. To get paragraphs numbered and counted, increase the default value of secnumdepth to 4

\titleformat {\section} [block] {\raggedright \fontsize{14}{14}\selectfont\bfseries} {\thesection.\space} {0pt} {}
\titlespacing {\section} {0pt} {12pt} {6pt}

\titleformat {\subsection} [block] {\raggedright \fontsize{12}{12}\selectfont\bfseries } {\thesubsection.\space} {0pt} {}
\titlespacing {\subsection} {0pt} {12pt} {6pt}

\titleformat {\subsubsection} [block] {\raggedright \fontsize{10}{10}\selectfont \bfseries} {\thesubsubsection.\space} {0pt} {}
\titlespacing {\subsubsection} {0pt} {12pt} {6pt}

\titleformat {\paragraph} [block] {\raggedright \fontsize{10}{10}\selectfont \bfseries} {} {0pt} {}
\titlespacing {\paragraph} {0pt} {12pt} {6pt}

%%%% Special section title style for back matter
\newcommand{\supplementary}[1]{
\par\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Supplementary Materials:} {#1}\par}}

\newcommand{\acknowledgments}[1]{
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Acknowledgments:} {#1}\par}}

\newcommand{\authorcontributions}[1]{%
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Author Contributions:} {#1}\par}}

\newcommand{\conflictsofinterest}[1]{%
\vspace{6pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Conflicts of Interest:} {#1}\par}}

\newcommand{\sampleavailability}[1]{%
\vspace{12pt}\noindent{\fontsize{9}{9}\selectfont\textbf{Sample Availability:} {#1}\par}}

\newcommand{\abbreviations}[1]{%
\vspace{12pt}\noindent{\selectfont\textbf{Abbreviations}\par\vspace{6pt}\noindent {\fontsize{9}{9}\selectfont #1}\par}}

%%%%% Defines the appendix
\renewcommand{\appendix}{%
\addcontentsline{toc}{section}{Appendix}
\setcounter{section}{0}%
\setcounter{subsection}{0}%
\setcounter{subsubsection}{0}%
%
\gdef\thesection{\@Alph\c@section}%
\gdef\thesubsection{\@Alph\c@section.\@arabic\c@subsection}%
%	

\gdef\theHsubsection{\@Alph\c@section.\@arabic\c@subsection}% for hyperref
\gdef\theHsubsection{\@Alph\c@section.\@arabic\c@subsection}% for hyperref
\csname appendixmore\endcsname
\renewcommand{\thefigure}{A\arabic{figure}}
\setcounter{figure}{0}
\renewcommand{\thetable}{A\arabic{table}}
\setcounter{table}{0}
\renewcommand{\thescheme}{A\arabic{scheme}}
\setcounter{scheme}{0}
\renewcommand{\thechart}{A\arabic{chart}}
\setcounter{chart}{0}
\renewcommand{\theboxenv}{A\arabic{boxenv}}
\setcounter{boxenv}{0}
\renewcommand{\theequation}{A\arabic{equation}}
\setcounter{equation}{0}
\renewcommand{\thetheorem}{A\arabic{theorem}}
\setcounter{theorem}{0}
\renewcommand{\thelemma}{A\arabic{lemma}}
\setcounter{lemma}{0}
\renewcommand{\thecorollary}{A\arabic{corollary}}
\setcounter{corollary}{0}
\renewcommand{\theproposition}{A\arabic{proposition}} 
\setcounter{proposition}{0} 
\renewcommand{\theexample}{A\arabic{example}}
\setcounter{example}{0} 
\renewcommand{\theremark}{A\arabic{remark}}
\setcounter{remark}{0} 
\renewcommand{\thedefinition}{A\arabic{definition}}
\setcounter{definition}{0} 
\renewcommand{\thehypothesis}{A\arabic{hypothesis}}
\setcounter{hypothesis}{0}
\renewcommand{\thenotation}{A\arabic{notation}}
\setcounter{notation}{0}
}


%%%% Layout
\RequirePackage[left=2.1cm,
				right=2.1cm,
				top=1.8cm,
				bottom=1.5cm,
				includehead,
				includefoot]{geometry}
\linespread{1.13}
\setlength{\parindent}{0.75cm}

%%%% Figures and tables
\newfloat{boxenv}{H}{lop}
\floatname{boxenv}{Box}

\newfloat{chart}{H}{lop}
\floatname{chart}{Chart}

\newfloat{scheme}{H}{lop}
\floatname{scheme}{Scheme}

\floatstyle{plain}
\restylefloat{figure}

\floatstyle{plaintop}
\restylefloat{table}

\floatstyle{plaintop}
\restylefloat{boxenv}

\floatstyle{plain}
\restylefloat{scheme}

\RequirePackage{caption} 
\captionsetup[figure]{position=bottom, margin=0cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

\captionsetup[scheme]{position=bottom, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

\captionsetup[chart]{position=bottom, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

\captionsetup[table]{position=top, margin=0cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

\captionsetup[boxenv]{position=top, margin=0.75cm, labelfont={bf, small, stretch=1.17}, labelsep=period, textfont={small, stretch=1.17}, aboveskip=6pt, belowskip=12pt, justification=justified}

%% For table footnotes
\newsavebox{\@justcentbox}
\newcommand{\justifyorcenter}[1]{
\sbox \@justcentbox{#1}
\ifdim \wd \@justcentbox >\hsize #1
\else \centerline{#1} \fi
}

%%%% Bullet lists
\newlength{\wideitemsep}
\setlength{\wideitemsep}{.5\itemsep}
\addtolength{\wideitemsep}{-7pt}
\let\olditem\item
\renewcommand{\item}{\setlength{\itemsep}{\wideitemsep}\olditem}

%%%% Header and footer (all pages except the first)
\renewcommand\headrule{} %% set line (from fancyhdr) in header to nothing
\newcommand{\changefont}{%
    \fontsize{8}{8}\selectfont
}
\pagestyle{fancy}
\fancyhead[L]{\changefont \MakeUppercase \leftmark }%section
\rhead{}
\lfoot{
\fontsize{8}{8} \selectfont {Version: {\@ \today}}
	}


%%%% Bibliography
\renewcommand\@biblabel[1]{#1.}
\def\thebibliography#1{
\linespread{1.44} 
\section*{References}
\addcontentsline{toc}{part}{References}
\fontsize{9}{9}\selectfont
\list{{\arabic{enumi}}}{\def\makelabel#1{\hss{#1}}
\topsep=0\p@
\parsep=5\p@
\partopsep=0\p@
\itemsep=0\p@
\labelsep=1.5mm
\settowidth\labelwidth{\footnotesize[#1]}%
\leftmargin\labelwidth
\advance\leftmargin\labelsep
\usecounter{enumi}}
}
\let\endthebibliography=\endlist

\endinput