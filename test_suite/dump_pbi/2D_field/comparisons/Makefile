PBI1_EXEC = ./old_ibmrun
PBI2_EXEC = ./ibmrun
SOURCECTRL = simpar1
CTRLFILE   = simpar
FIELD = $(shell grep output_var simpar1 | tr "[[:lower:]]" "[[:upper:]]"| sed 's/OUTPUT_VAR//;;s/!.*//;;s/=//;;s/ //g' )

.PHONY = diff one two variables clean all
diff:
	make one two
	compare $(FIELD)_one.png $(FIELD)_two.png $(FIELD)_image_diff.png
	R --slave < Plot_2D_field_diff.r
	convert field_diff.pdf $(FIELD)_field_diff.png 

one: $(CTRLFILE) 
	$(PBI1_EXEC) $(CTRLFILE)
	R --slave < Plot_2D_field.r
	convert $(FIELD)_2D_fields.pdf $(FIELD)_one.png
	mv 2D_fields.txt one_fields.txt 

two: $(CTRLFILE) 
	$(PBI2_EXEC) $(CTRLFILE)
	R --slave < Plot_2D_field.r
	convert $(FIELD)_2D_fields.pdf $(FIELD)_two.png
	mv 2D_fields.txt two_fields.txt 

$(CTRLFILE): $(SOURCECTRL)
	cp $(SOURCECTRL) $(CTRLFILE)

variables:
	@echo "FIELD = $(FIELD)"

clean:
	-rm *.png *.pdf *.txt
	-rm $(CTRLFILE)

all:
	make clean
	make depth temperature salinity U V hderiv vderiv hdiffus vdiffus 
	convert *.png comparison.pdf

%: 
	sed 's/^output_var.*/output_var   =$*/' $(SOURCECTRL) > $(CTRLFILE)
	make diff
	
