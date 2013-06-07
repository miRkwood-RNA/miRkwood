sub input_file {
    return File::Spec->catfile($FindBin::Bin, 'data', @_);
}

1;
