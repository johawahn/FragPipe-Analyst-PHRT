$(document).on("click", ".clickable-image", function() {
  var imgSrc = $(this).find("img").attr("src");
  Shiny.setInputValue("click_image", imgSrc);
});

