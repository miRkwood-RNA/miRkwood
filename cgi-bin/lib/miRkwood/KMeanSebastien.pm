# ABSTRACT: Handling the cluster making from a BAM file.

package KMeanSebastien;

use strict;
use warnings;

use constant {
	FIRST_CLASS_CREATED => 1,
	SECOND_CLASS_CREATED => 2,
	ASSIGNED_ONLY_EXISTING_CLASS => 3,
	ASSIGNED_FIRST_CLASS => 4,
	ASSIGNED_SECOND_CLASS => 5
};

=method new

Constructor

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {
		class_1_total => 0,
		class_1_count => 0,
		class_1_avg => 0,
		class_0_total => 0,
		class_0_count => 0,
		class_0_avg => 0,
		class_border => 0,
		has_different_values => 0,
		points => []
    }, $class;
    return $self;
}

sub __binary_insert {
    my ($ary, $val) = @_;
    my ($min, $max, $middle, $last) = (0, scalar @{$ary}, 0, -1);
    if ($max > 0) {
		while (1) {
			$middle = int(($min+$max)/2);
			if ($min == $max) {
				last;
			}
			if ($ary->[$middle] < $val) {
				$min = $middle+1;
			}
			elsif ($ary->[$middle] > $val) {
				$max = $middle;
			}
			else {
				last;
			}
		}
    }
    splice @$ary, $middle, 0, $val;
    return $middle;
}


sub __recompute_avg {
	my $this = shift;

	$this->{class_0_avg} = $this->{class_0_total}/$this->{class_0_count};
	$this->{class_1_avg} = $this->{class_1_total}/$this->{class_1_count};

	while (1) {
		# Assignment step
		# In this step, for every point we check wether it is closer to the class 0 avg or the class 1 avg, and we reassign it if needed.
		# We try to do it smartly. As we are trying to find two clusters in 1D, the two clusters are speratated by only one point/margin.
		# This means that only the points near the margin are likely to be reassigned.
		# We start to check the points near the margin upstream. If nothing gets changed, we check downstream. Note that it's impossible that we must reassign on both directions.
		# That's because the margin will either move left or right but not both.

		# Check near the border, in the class 1 cluster (i.e. try to move the margin right)
		my $assigned_upstream = 0;
		for (my $i = $this->{class_border}, my $e = scalar(@{$this->{points}}); $i < $e; $i++) {
			if (abs($this->{points}[$i] - $this->{class_0_avg}) < abs($this->{points}[$i] - $this->{class_1_avg})) {
				$assigned_upstream++;
				$this->{class_0_total} += $this->{points}[$i];
				$this->{class_1_total} -= $this->{points}[$i];
# 				$this->{class_0_count}++;
# 				$this->{class_1_count}--;

# 				$this->{class_0_avg} = $this->{class_0_total}/$this->{class_0_count};
# 				$this->{class_1_avg} = $this->{class_1_total}/$this->{class_1_count};
			}
			else {
				last;
			}
		}
		# If we couldn't move the margin right, we try to move it left
		if ($assigned_upstream == 0) {
			for (my $i = $this->{class_border}-1; $i >= 0; $i--) {
				if (abs($this->{points}[$i] - $this->{class_0_avg}) > abs($this->{points}[$i] - $this->{class_1_avg})) {
					$assigned_upstream--;
					$this->{class_0_total} -= $this->{points}[$i];
					$this->{class_1_total} += $this->{points}[$i];
					
# 					$this->{class_0_count}--;
# 					$this->{class_1_count}++;
# 
# 					$this->{class_0_avg} = $this->{class_0_total}/$this->{class_0_count};
# 					$this->{class_1_avg} = $this->{class_1_total}/$this->{class_1_count};
				}
				else {
					last;
				}
			}
		}
		# If nothing has changed (nothing was misclassified), the algorithm converged
		if ($assigned_upstream == 0) {
			last;
		}
		$this->{class_0_count} += $assigned_upstream;
		$this->{class_1_count} -= $assigned_upstream;
		# Recompute the averages
		$this->{class_0_avg} = $this->{class_0_total}/$this->{class_0_count};
		$this->{class_1_avg} = $this->{class_1_total}/$this->{class_1_count};
		$this->{class_border} = $this->{class_0_count};
    }
}


# Returns a hash (chr_name => read positions)
sub add_point {
    my ($this, $value) = @_;
    my $old_val = $this->{points}[0];
	my $i = __binary_insert($this->{points}, $value);
	if (scalar @{$this->{points}} == 1) {
		$this->{class_0_count} = 1;
		$this->{class_0_avg} = $value;
		$this->{class_0_total} = $value;
		$this->{class_border} = 1;
		return FIRST_CLASS_CREATED;
	}
	elsif ($this->{has_different_values} == 0) {
		if ($value != $old_val) {
			$this->{has_different_values} = 1;
			$this->{class_1_count} = 1;
			$this->{class_1_total} = $value;
			$this->{class_1_avg} = $value;
			$this->{class_border} = $i;
			if ($value < $old_val) {
				my $tmp = $this->{class_0_count};
				$this->{class_0_count} = $this->{class_1_count};
				$this->{class_1_count} = $tmp;

				$tmp = $this->{class_0_total};
				$this->{class_0_total} = $this->{class_1_total};
				$this->{class_1_total} = $tmp;

				$tmp = $this->{class_0_avg};
				$this->{class_0_avg} = $this->{class_1_avg};
				$this->{class_1_avg} = $tmp;

				$this->{class_border} = 1;
			}
			return SECOND_CLASS_CREATED;
		}
		else {
			$this->{class_0_count}++;
			$this->{class_0_total} += $value;
			$this->{class_border}++;
			return ASSIGNED_ONLY_EXISTING_CLASS;
		}
	}
	else {
		if (abs($value - $this->{class_0_avg}) < abs($value - $this->{class_1_avg})) {
			$this->{class_border}++;
			$this->{class_0_count}++;
			$this->{class_0_total} += $value;
		}
		else {
			$this->{class_1_count}++;
			$this->{class_1_total} += $value;
		}
		# Recompute the avg
		$this->__recompute_avg();
		return $i < $this->{class_border} ? ASSIGNED_FIRST_CLASS : ASSIGNED_SECOND_CLASS;
	}
}

# value should be in the set
sub class_of {
	my ($this, $value) = @_;
	if ($this->{class_1_count} > 0) {
		return $value < $this->{points}[$this->{class_border}] ? ASSIGNED_FIRST_CLASS : ASSIGNED_SECOND_CLASS;
	}
	return ASSIGNED_ONLY_EXISTING_CLASS;
}

sub clear_points {
	my $this = shift;
	$this->{class_1_total} = 0;
	$this->{class_1_count} = 0;
	$this->{class_0_count} = 0;
	$this->{class_0_total} = 0;
	$this->{class_border} = 0;
	$this->{has_different_values} = 0;
	$this->{points} = [];
}

1;
