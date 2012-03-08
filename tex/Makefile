LATEX       = pdflatex
RM          = rm -rf
RM_TMP      = ${RM} $(foreach suff, ${TMP_SUFFS}, ${NAME}.${suff})

TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out
SUFF        = pdf

CHECK_RERUN = grep Rerun $*.log

NAME    = transits
DOC_OUT = ${NAME}.${SUFF}

default: ${DOC_OUT}

%.pdf: %.tex
	${LATEX} $<
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."

plots:
	$(foreach f, plotting/*.py, python ${f};)

clean:
	${RM_TMP}

