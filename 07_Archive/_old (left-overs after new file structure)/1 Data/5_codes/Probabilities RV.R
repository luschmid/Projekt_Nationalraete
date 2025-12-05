

part_unelected=0.6
part_elected=0.75

el_unelected_cond=0.5
el_elected_cond=0.76

el_unelected_uncond=0.27
el_elected_uncond=0.6

el_unelected_cond*0.12/(el_unelected_cond*0.12+part_unelected*0.27)

log(part_elected)-log(part_unelected)+log(el_elected_cond-log(el_unelected_cond))
log(el_elected_uncond)-log(el_unelected_uncond)


part_unelected=0.6
part_elected=0.75
el_unelected_cond=0.5
el_elected_cond=0.5
el_unelected_uncond=0.3
el_elected_uncond=0.375