def myvariant_wrapper(query_term):
    import myvariant
    mv = myvariant.MyVariantInfo()
    return mv.getvariant(query_term)
