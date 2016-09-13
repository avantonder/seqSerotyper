# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: tfasta.pm,v 1.16 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA3::tfasta;

use NPB::Parse::Format::FASTA3;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::HEADER);


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA::RANK);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    #ranked search hits
    while (defined ($line = $text->next_line)) {
	
	next    if $line =~ /$NPB::Parse::Format::FASTA3::RANK_START/o;

	#pre-tfast[axy]3.4 behaviour
	if($line =~ /^
	   \s*
	   (\S+)                #id
	   \s+
	   (.*)                 #desc (may be empty)
	   \s+
	   \(\s*(\d+)\)         #aa
	   \s+
	   \[(\S)\]             #frame
	   \s+
	   (\d+)                #initn
	   \s+
	   (\S+)                #init1
	   \s+
	   (\d+)                #opt
	   \s+
	   (\S+)                #z-score
	   \s+
	   (\S+)                #E(205044)
	   \s*
	   $/xo) {
	    
	    $self->test_args($line, $1, $3,$4, $5,$6,$7,$8,$9); #not $2
	    
	    push(@{$self->{'hit'}},
		 { 
		  'id'     => NPB::Parse::Record::clean_identifier($1),
		  'desc'   => $2,
		  'length' => $3,
		  'frame'  => NPB::Parse::Format::FASTA::parse_frame($4),
		  'orient' => NPB::Parse::Format::FASTA::parse_orient($4),
		  'initn'  => $5,
		  'init1'  => $6,
		  'opt'    => $7,
		  'zscore' => $8,
		  'expect' => $9,
		 });
	    next;
	}
    
	#tfast[axy]3.4
	if($line =~ /^
	   \s*
	   (\S+)                #id
	   \s+
	   (.*)                 #desc (may be empty)
	   \s+
	   \(\s*(\d+)\)         #aa
	   \s+
	   \[(\S)\]             #frame
	   \s+
	   (\d+)                #headed 'initn' but actually 'opt'
	   \s+
	   (\S+)                #headed 'init1' but actually 'bits'
	   \s+
	   (\S+)                #E(205044)
	   \s*
	   $/xo) {
	    
	    $self->test_args($line, $1, $3,$4, $5,$6,$7); #not $2
	    
	    push(@{$self->{'hit'}},
		 { 
		  'id'     => NPB::Parse::Record::clean_identifier($1),
		  'desc'   => $2,
		  'length' => $3,
		  'frame'  => NPB::Parse::Format::FASTA::parse_frame($4),
		  'orient' => NPB::Parse::Format::FASTA::parse_orient($4),
		  'initn'  => '',
		  'init1'  => '',
		  'opt'    => $5,
		  'bits'   => $6,
		  'zscore' => '',
		  'expect' => $7,
		 });
	    next;
	}
    
	#blank line or empty record: ignore
	next    if $line =~ /$NPB::Parse::Format::FASTA3::NULL/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::TRAILER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::MATCH::SUM;

use vars qw(@ISA);
use NPB::Parse::Regexps;

@ISA = qw(NPB::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    while (defined ($line = $text->next_line(1))) {

	next  if $line =~ /^\s*$/;

	if ($line =~ /^>+/) {
	    $record = $line;       #process this later
	    next;
	}

	#3.1t07; known to work with tfastx
	if ($line =~ /^
	    (rev-comp)?            #frame
	    \s*
	    initn\:\s*(\S+)        #initn
	    \s*
	    init1\:\s*(\S+)        #init1
	    \s*
	    opt\:\s*(\S+)          #opt
	    \s*
	    Z-score\:\s*(\S+)      #z
	    \s*
	    expect\(\)\s*(\S+)     #E
	    \s*
	    $/xo) {
	    
	    $self->test_args($line,$2,$3,$4,$5,$6);
	    
	    (
	     $self->{'frame'},
	     $self->{'orient'},
	     $self->{'initn'},
	     $self->{'init1'},
	     $self->{'opt'},
	     $self->{'zscore'},
	     $self->{'bits'},
	     $self->{'expect'},
	    ) = (
		NPB::Parse::Format::FASTA::parse_frame($1),
		NPB::Parse::Format::FASTA::parse_orient($1),
		$2, $3, $4, $5, '', $6,
	    );
	    next;
	}

	#3.4t23; known to work with tfasta
	if ($line =~ /^
	    Frame\:\s*(\S+)        #frame 1-6
	    \s+
	    initn\:\s*(\S+)        #initn
	    \s+
	    init1\:\s*(\S+)        #init1
	    \s+
	    opt\:\s*(\S+)          #opt
	    \s+
	    Z-score\:\s*(\S+)      #z
	    \s+
	    bits\:\s*(\S+)         #bits
	    \s+
	    E\(\)\:\s*(\S+)        #E
	    \s*
	    $/xo) {

	    $self->test_args($line,$1,$2,$3,$4,$5,$6,$7);

	    (
	     $self->{'frame'},
	     $self->{'orient'},
	     $self->{'initn'},
	     $self->{'init1'},
	     $self->{'opt'},
	     $self->{'zscore'},
	     $self->{'bits'},
	     $self->{'expect'},
	    ) = (
		NPB::Parse::Format::FASTA::parse_frame($1),
		NPB::Parse::Format::FASTA::parse_orient($1),
		$2, $3, $4, $5, $6, $7,
	    );
	    next;
	}

	#3.1t07; known to work with tfastx
	if ($line =~ /^
	    (?:Smith-Waterman\s+score:\s*(\d+);)?    #sw score
	    \s*($RX_Ureal)%                          #percent identity
	    \s*identity\s+in\s+(\d+)                 #overlap length
	    \s+(?:aa|nt)\s+overlap
	    (?:\s+\((\S+)\))?                        #sequence ranges
	    /xo) {

	    $self->test_args($line,$2,$3);
	    
	    (
	     $self->{'score'},
	     $self->{'id_percent'},
	     $self->{'ungap_percent'},
	     $self->{'overlap'},
	     $self->{'ranges'},
	    ) = (defined $1?$1:0,$2,'',$3,defined $4?$4:'');
	    next;
        }

	#3.4t23; known to work with tfasta
	if ($line =~ /^
	    (?:(?:banded\s+)?Smith-Waterman\s+score:\s*(\S+);)?  #sw score
	    \s*(\S+)%\s+identity                     #percent identity
	    \s+\((\S+)%\s+ungapped\)                 #percent ungapped
	    \s+in\s+(\d+)                            #overlap length
	    \s+(?:aa|nt)\s+overlap
	    \s*
	    (?:\s+\((\S+)\))?                        #sequence ranges
	    /xo) {

	    $self->test_args($line,$2,$3);
	    
	    (
	     $self->{'score'},
	     $self->{'id_percent'},
	     $self->{'ungap_percent'},
	     $self->{'overlap'},
	     $self->{'ranges'},
	    ) = (defined $1?$1:0,$2,$3,$4,defined $5?$5:'');
	    next;
        }

	#should only get here for multiline descriptions
	if (defined $self->{'initn'}) {
	    $self->warn("unknown field: $line");
	    next;
	}

	#accumulate multiline descriptions (fasta... -L)
	$record .= ' ' . $line;
    }

    #now split out the description
    if ($record =~ /^
	>+
	(\S+)                      #id
	\s+
	(.*)                       #description
	\s+
	\(\s*(\d+)\s*(?:aa|nt)\)   #length
	\s*
	$/xo) {

	$self->test_args($record, $1, $2, $3);

	(
	 $self->{'id'},
	 $self->{'desc'},
	 $self->{'length'},
	) = (NPB::Parse::Record::clean_identifier($1),
	     NPB::Parse::Record::strip_english_newlines($2), $3);
    } else {
	$self->warn("unknown field: $record");
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3::tfasta::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::MATCH::ALN);

# tfast[axy]  pro x dna
#
# program		sbjct_base	range in summary?
# tfastx  3.0t82	1		no
# tfastxy 3.1t07	1		no
# tfasta  3.4t23	3		yes
# tfastx  3.4t23	3		yes
# tfasty  3.4t23	3		yes

sub query_orient { '+' }
sub sbjct_orient { $_[0]->get_summary->{'orient'} }

sub query_base { return 1 }
sub sbjct_base {
    return 1  if ($_[0]->get_header->{'version'} cmp "3.4") < 0;
    return 3;
}


###########################################################################
1;
