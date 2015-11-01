SRC=$(wildcard *.md) 
OBJS=$(SRC:.md=.html)

slides: $(OBJS)

$(OBJS) : %.html : %.md 
	#pandoc -s --mathjax -t slidy --highlight-style pygments session1.md -o session1.html
	pandoc -s --mathjax -t slidy --highlight-style pygments $< -o $@



