;   14/03/2018  JARM
;
;   Given an input flag table, returns a mask where a certain bit is active
;

PRO sancho_check_activebit_inflags, flag, bit, masc
if n_params() lt 2 then begin
	print,""
	print,"  Syntax -- sancho_check_activebit_inflags, flag, bit, masc "
	print,""
	print,"  INPUTS: "
	print,"      flag(nchan,ndata) table"
	print,"      bit number"
	print,""
	print,"  OUTPUTS:"
	print,"      masc(nchan,ndata) showing where the bit is active."
	print,""
	return
endif
N    = bit-1L
masc = flag/2^N mod 2
END
