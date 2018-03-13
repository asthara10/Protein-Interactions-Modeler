#!/bin/perl
#

 use Getopt::Long;

 
 
 &GetOptions("i=s" => \$input,
             "d=s" => \@database,
             "h" => \$help,
             "o=s" => \$output );


   if (defined $help){
      printf " 
               i input of codes
               d database of sequences
               h print help
               o output 
                \n";
   }

     open (KIN,"<$input");
     while(<KIN>){
      @query=split;
      $name={};
      $name={
	      swall=>$query[0],
	      all  =>[@query]};
      push @code,$name;
      
     }
     close KIN;

     for $j (0..$#database){

     open (SWP,"cat $database[$j]|");
     $name1=" ";
     $name2=" ";
     undef @seq;
     while(<SWP>){
      if (/^>(\S+)/){
	      $query{$name1}=join "",@seq;
	      $query{$name2}=join "",@seq;
	      $name1=$1;
	      $name2=$1;
              if (/^>\S+ \((\S+)\) /){$name2=$1;}
	      undef @seq;
	      next;}
      chop;

      push @seq,$_;
     }
     close SWP;

     }


     open OUT,">$output";
     for $prot (@code){
          @residue=split ("",$query{$prot->{swall}});
          if ($#residue>0) {printf OUT ">%s %s\n",$prot->{all}->[0],$prot->{all}->[1];}
          $n=0;
          for $res (@residue){
              if ($n>59){print OUT "\n";$n=0}
              print OUT "$res";
              $n++;
          }
          if ($#residue>0) {print OUT "\n";}
     }
 
