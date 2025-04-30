% f_activator
function f=Blmp1_logic(H_dafREblmp,H_lin42ACblmp,H_lin29REblmp,H_blmpACblmp,H_dre1ACblmp,H_dafACblmp,setNo)
 
if setNo == 2
    f = ((H_lin42ACblmp *H_dafREblmp *H_lin29REblmp) +H_blmpACblmp) -((H_lin42ACblmp *H_dafREblmp *H_lin29REblmp) * H_blmpACblmp);
	%f = (H_lin42ACblmp *H_dafREblmp *H_lin29REblmp);
	%f = (H_lin42ACblmp *H_dafREblmp *H_lin29REblmp*H_blmpACblmp)

elseif setNo == 1
    f = ((H_lin42ACblmp *H_dafREblmp *H_lin29REblmp) * H_blmpACblmp);
   %f = (H_lin42ACblmp *H_dafREblmp *H_lin29REblmp);
elseif setNo == 7
    f = ((((H_lin42ACblmp + H_dafREblmp)- (H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp) +H_blmpACblmp) -((((H_lin42ACblmp +H_dafREblmp)- (H_lin42ACblmp * H_dafREblmp)) *H_lin29REblmp) *H_blmpACblmp);
elseif setNo == 4
    f = ((((H_lin42ACblmp +H_dafREblmp)- (H_lin42ACblmp * H_dafREblmp)) *H_lin29REblmp) *H_blmpACblmp);
elseif setNo == 8
    f = (((((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) +H_blmpACblmp) -(((((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) *H_blmpACblmp); 
elseif setNo == 5
    f = (((((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp +H_dafREblmp) -(H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) *H_blmpACblmp); 
elseif setNo == 6
    f = (((((H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) +H_blmpACblmp) -(((((H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) *H_blmpACblmp); 
elseif setNo == 3
    f = (((((H_lin42ACblmp *H_dafREblmp)) +H_lin29REblmp)-(((H_lin42ACblmp *H_dafREblmp)) *H_lin29REblmp)) *H_blmpACblmp); 
end
end
