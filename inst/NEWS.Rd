% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/diptest/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{diptest}}
\encoding{UTF-8}

\section{CHANGES in diptest VERSION 0.76-0 (2021-03-23)}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item add \file{README.Rd} mostly for github readers
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Added Imports (to "base" packages where they were not checked previously) to \file{NAMESPACE}.
      \item Using \file{NEWS.Rd} file more.
    }
  }
}

\section{CHANGES in diptest VERSION 0.75-7 (2015-06-07)}{% CRAN release

  \subsection{NEW FEATURES}{
    \itemize{
      \item Started this \file{NEWS.Rd} file, to eventually replace the \file{ChangeLog}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{CHANGES in diptest VERSION 0.75-6 (2014-11-25)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{CHANGES in diptest VERSION 0.75-5 (2013-07-23)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item add \file{NEWS.Rd} (albeit mostly empty)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{rdRDS()} wrapper corrrectly tests for R 2.13.0
    }
  }
}


%% this is the *latest* entry  in  ../ChangeLog
\section{CHANGES in diptest VERSION 0.75-4 (2012-08-13)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Enable package for pre-R-2.13.x via \code{rdRDS()} wrapper

      \item \code{dip.test()} now also returns an \code{alternative} component,
      e.g. for printing.
    }
  }
}

\section{CHANGES in diptest VERSION 0.75-3 (2012-04-18)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


\section{CHANGES in diptest VERSION 0.75-1 (2011-08-10)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


\section{CHANGES in diptest VERSION 0.25-3 (2010-08-11)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item First version of the \dQuote{"diptest issues"} vignette.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


\section{CHANGES in diptest VERSION 0.25-2 (2009-02-09)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


\section{CHANGES in diptest VERSION 0.25-1 (2004-08-12)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


\section{CHANGES in diptest VERSION 0.25-0 (2004-02-13)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item More output in the \R object, allows \code{debug} information.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Previously, the \code{dip()} had not been \dQuote{symmetric}
      with respect to \dQuote{mirroring} such as a sign flip in the
      data.  Thanks to Yong Lu, who was able to track the bug to a
      misplaced \code{(} in the original Fortran code, the bug has been
      fixed on Oct.10, 2003.
    }
  }
}

\section{CHANGES in diptest VERSION 0.9-1 (2003-07-15)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item First(?) CRAN release
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item on 1994-07-30, added code to prevent an infinite loop in
      rare cases, but e.g., for \code{dip(1:n)}.
    }
  }
}

