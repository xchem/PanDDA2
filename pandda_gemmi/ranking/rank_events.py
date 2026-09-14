def rank_events(events, sites, autobuilds, rank_method, existing_sites):

    ranked_event_ids = rank_method(events, sites, autobuilds, existing_sites)

    return ranked_event_ids
