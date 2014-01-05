#
                                                   sub j(\$){($
                     P,$V)=                      @_;while($$P=~s:^
                 ([()])::x){                    $V+=('('eq$1)?-32:31
           }$V+=ord(  substr(                 $$P,0,1,""))-74} sub a{
          my($I,$K,$  J,$L)=@_               ;$I=int($I*$M/$Z);$K=int(
         $K*$M/$Z);$J=int($J*$M             /$Z);$L=int($L*$M/$Z); $G=$
         J-$I;$F=$L-$K;$E=(abs($          G)>=abs($F))?$G:$F;($E<0) and($
          I,$K)=($J,$L);$E||=.01       ;for($i=0;$i<=abs$E;$i++ ){ $D->{$K
                  +int($i*$F/$E)      }->{$I+int($i*$G/$E)}=1}}sub p{$D={};$
                 Z=$z||.01;map{    $H=$_;$I=$N=j$H;$K=$O=j$H;while($H){$q=ord
                substr($H,0,1,"" );if(42==$q){$J=j$H;$L=j$H}else{$q-=43;$L =$q
              %9;$J=($q-$L)/9;$L=$q-9*$J-4;$J-=4}$J+=$I;$L+=$K;a($I,$K,$J,$ L);
              ($I,$K)=($J,$L)}a($I,$K,$N,$O)}@_;my$T;map{$y=$_;map{ $T.=$D->{$y}
              ->{$_}?$\:' '}(-59..59);$T.="\n"}(-23..23);print"\e[H$T"}$w= eval{
              require Win32::Console::ANSI};$b=$w?'1;7;':"";($j,$u,$s,$t,$a,$n,$o
              ,$h,$c,$k,$p,$e,$r,$l,$C)=split/}/,'Tw*JSK8IAg*PJ[*J@wR}*JR]*QJ[*J'.
               'BA*JQK8I*JC}KUz]BAIJT]*QJ[R?-R[e]\RI'.'}Tn*JQ]wRAI*JDnR8QAU}wT8KT'.
               ']n*JEI*EJR*QJ]*JR*DJ@IQ[}*JSe*JD[n]*JPe*'.'JBI/KI}T8@?PcdnfgVCBRcP'.
                '?ABKV]]}*JWe*JD[n]*JPe*JC?8B*JE};Vq*OJQ/IP['.'wQ}*JWeOe{n*EERk8;'.
                  'J*JC}/U*OJd[OI@*BJ*JXn*J>w]U}CWq*OJc8KJ?O[e]U/T*QJP?}*JSe*JCnTe'.
                   'QIAKJR}*JV]wRAI*J?}T]*RJcJI[\]3;U]Uq*PM[wV]W]WCT*DM*SJ'.  'ZP[Z'.
                      'PZa[\]UKVgogK9K*QJ[\]n[RI@*EH@IddR[Q[]T]T]T3o[dk*JE'.  '[Z\U'.
                        '{T]*JPKTKK]*OJ[QIO[PIQIO[[gUKU\k*JE+J+J5R5AI*EJ00'.  'BCB*'.
                             'DMKKJIR[Q+*EJ0*EK';sub h{$\ = qw(% & @ x)[int    rand
                              4];map{printf  "\e[$b;%dm",int(rand 6)+101-60*   ($w
                               ||0);system(  "cls")if$w ;($A,$S)=    ($_[1],   $
                                _[0]);($M,   @,)= split  '}';for(     $z=256
                                ;$z>0; $z   -=$S){$S*=   $A;p @,}      sleep$_
                                [2];while   ($_[3]&&($    z+=$ S)       <=256){
                                p@,}}("".   "32}7D$j"     ."}AG".       "$u}OG"
                                ."$s}WG"    ."$t",""      ."24}("        ."IJ$a"
                                ."}1G$n"    ."}CO$o"     ."}GG$t"        ."}QC"
                                 ."$h}"      ."^G$e"    ."})IG"          ."$r",
                                 "32}?"       ."H$p}FG$e}QG$r".          "}ZC"
                                 ."$l",          "28}(LC" .""            ."".
                                 "$h}:"           ."J$a}EG".             "$c"
                                 ."}M"             ."C$k}ZG".            "$e"
                                 ."}"             ."dG$r","18"          ."}("
                                ."D;"            ."$C"  )}{h(16         ,1,1,0
                               );h(8,          .98,0,0   );h(16         ,1,1,1)
                               ;h(8.0         ,0.98,0,     1);         redo}###
                             #written                                 060204 by
                           #liverpole                                  @@@@@@@
                        #@@@@@@@@@@@

