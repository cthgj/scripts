# DEFINE VARIABLES
#-----------------
annot_T=$1
onto1=$2
onto2=$3
error="Usage: $0 <annotation table> <ontology 1> [ ontology 2 ]"
dbname="obo";


if ([ -z $annot_T ] || [ -z $onto1 ]) then 
echo $error
exit
fi


if ([ -z $onto2 ]) then
onto2=$onto1;
fi;

o1=`echo $onto1 | tr "." "_"`;
o2=`echo $onto2 | tr "." "_"`;
table="temp_"$annot_T"_"$o1
echo $table

# DROP TEMPORARY TABLES
#----------------------
# echo "DROP TABLE temp;"             | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
# echo "DROP TABLE temp_1;"           | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_asso;"        | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_count_asso;"  | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_count_goid;"  | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_nb_genes;"    | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_good_genes;"  | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE temp_final;"       | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;
 echo "DROP TABLE $table;"           | psql  -U cchapple -h 10.1.1.53 -q -d $dbname;

# FETCH DIRECT ANNOTATIONS
#-------------------------
#if [ "$onto1" == "$onto2" ] 
#then
# QUERY="select T1.dbid, T1.goid as goid1, T2.goid as goid2 into temp table temp from $annot_T T1, $annot_T T2 where T1.dbid=T2.dbid AND T1.goid < T2.goid AND T1.subonto='$onto1' AND T2.subonto='$onto2' group by T1.dbid, goid1, goid2;"
#else
# QUERY="select T1.dbid, T1.goid as goid1, T2.goid as goid2 into temp from $annot_T T1, $annot_T T2 where T1.dbid=T2.dbid  AND T1.subonto='$onto1' AND T2.subonto='$onto2' group by T1.dbid, goid1, goid2;"
#fi



#QUERY="CREATE TABLE temp_asso (gene char(50), id1 char(30), id2 char(30));"

#echo $QUERY | psql -q -d $dbname 

## This will sort the table in order to avoid A<=>B and B<=>A
## each gene contains all the combinations of its direct annotations
## and the implicit ancestor annotations
## This is used to find the no of genes annotated to any pair of GOs.

QUERY="select T1.dbid as gene, (case when ( T4.ancestor_id < T5.ancestor_id) then T4.ancestor_id else T5.ancestor_id end) as id1, (case when ( T4.ancestor_id >= T5.ancestor_id) then T4.ancestor_id else T5.ancestor_id end) as id2 into table temp_asso from $annot_T T1, $annot_T T2, ancestors T4, ancestors T5 where T1.dbid=T2.dbid AND T4.id=T1.goid AND T5.id=T2.goid AND T1.subonto='$onto1' AND T2.subonto='$onto2' group by gene,id1,id2;"

# QUERY="select T1.dbid as gene, 
# (case when ( T4.ancestor_id < T5.ancestor_id) then 
# T4.ancestor_id else T5.ancestor_id end) as id1, 
# (case when ( T4.ancestor_id >= T5.ancestor_id) then 
# T4.ancestor_id else T5.ancestor_id end) as id2
# from $annot_T T1, $annot_T T2, ancestors T4, ancestors T5 
# where T1.dbid=T2.dbid AND T4.id=T1.goid AND T5.id=T2.goid AND T1.subonto='$onto1' AND T2.subonto='$onto2' group by gene,id1,id2;"
echo "******************000******************************"
echo $QUERY
echo $QUERY | psql  -U cchapple -h 10.1.1.53 -q -d $dbname 




# FETCH ANCESTORS
#----------------

# echo "Now fetching ancestors..."

# QUERY="select T1.dbid, T3.ancestor_id, T1.goid2 into temp_1 from ancestors T3, 
# temp T1 where T1.goid1=T3.id  
# AND T3.is_alt_id=FALSE ; "

# echo $QUERY
# echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname  

# QUERY="select T1.dbid, T1.ancestor_id as a1, T2.ancestor_id as a2 into temp_asso from temp_1 T1, ancestors T2 where T1.goid2=T2.id AND T2.is_alt_id=FALSE group by T1.dbid, a1, a2"

# echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

QUERY="select count(gene) as nb_asso, id1 as a1, id2 as a2 into table temp_count_asso from temp_asso group by a1, a2;"
echo "******************111******************************"

echo $QUERY
echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 


# count how many genes have annotations in two distinct ontologies: (this is the reference set)
#==================================================================

if [ "$onto1" != "$onto2" ] 
then

    QUERY="select count(distinct T1.dbid) as nb_good  into table temp_nb_genes from  $annot_T T1, (select distinct dbid from $annot_T where subonto='$onto1') T2
where T2.dbid=T1.dbid AND T1.subonto='$onto2';"

    echo $QUERY
    echo "******************222******************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

    QUERY="select distinct T1.dbid as good  into table temp_good_genes from  $annot_T T1, (select distinct dbid from $annot_T where subonto='$onto1') T2
where T2.dbid=T1.dbid AND T1.subonto='$onto2';"

    echo $QUERY
    echo "*****************333*******************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

else

# count how many genes have >=TWO annotations in one ontology:
#=============================================================

    QUERY="select count(T1.dbid) as nb_good into table temp_nb_genes from 
(select dbid, count(distinct goid) as nb from $annot_T where subonto='$onto1' group by dbid) T1 where nb >1;"

    echo $QUERY
    echo "******************444******************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

    QUERY="select T1.dbid as good into table temp_good_genes from 
(select dbid, count(distinct goid) as nb from $annot_T where subonto='$onto1' group by dbid) T1 where nb >1;"

    echo $QUERY
    echo "******************555******************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

fi


# NOW COUNT OCCURENCES OF SINGLE GOID
#------------------------------------

if [ "$onto1" == "$onto2" ] 
then

    QUERY="select count(distinct T1.dbid) as nb_goid, T2.ancestor_id into table temp_count_goid from $annot_T T1, ancestors T2, temp_good_genes T3 where T1.subonto='$onto1' AND T1.goid = T2.id AND T1.dbid=T3.good group by T2.ancestor_id ;"
    
    echo $QUERY
    echo "******************666*****************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 
else
    QUERY="select count(distinct T1.dbid) as nb_goid, T2.ancestor_id into table temp_count_goid from $annot_T T1, ancestors T2, temp_good_genes T3 where (T1.subonto='$onto2' OR T1.subonto='$onto1') AND T1.goid = T2.id AND T1.dbid=T3.good group by T2.ancestor_id ;"
    
    echo $QUERY
    echo "*******************777*****************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 
fi



echo "table="$table

# NOW OUTPUT RESULTS
#-------------------

QUERY="SELECT T1.a1 as id1,
	T1.a2 as id2, 
        common_anc(T1.a1,T1.a2),
	T2.nb_good as nb_genes, 
	T3.nb_goid as nb_id1, 
	T4.nb_goid as nb_id2, 
	T1.nb_asso as nb_both, 
	T5.term as term1, 
	T6.term as term2 into table temp_final from 
			temp_count_asso T1, 
			temp_nb_genes T2, 
			temp_count_goid T3, 
			temp_count_goid T4, 
			terms T5, terms T6
	where T1.a1=T3.ancestor_id 
	AND T1.a2 = T4.ancestor_id 
	AND T1.a1=T5.id 
	AND T1.a2=T6.id 
		order by nb_both desc;"

echo $QUERY
echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname 

# the last query is meant to eliminate combinations (id1,id2) where one is the ancestor of the other
# we use a EXCEPT clause to remove these combinations.

echo "Remove direct filiation relationships..."

QUERY="SELECT id1, id2, common_anc, nb_genes, nb_id1,nb_id2,nb_both, term1, term2 into $table from temp_final
except select T1.id1,T1.id2,T1.common_anc,T1.nb_genes,T1.nb_id1,T1.nb_id2,T1.nb_both, T1.term1, T1.term2 from temp_final T1, 
ancestors T2
where ( (T1.id1 = T2.id AND T1.id2=T2.ancestor_id) OR 
(T1.id2 = T2.id AND T1.id1=T2.ancestor_id));"

echo $QUERY;
    echo "*****************888*******************************"
    echo $QUERY | psql -U cchapple -h 10.1.1.53  -q -d $dbname;


echo "Done!"



#### nb_asso is a list of all GOs and the number of genes annotated to each of them (implicit and explicit)
